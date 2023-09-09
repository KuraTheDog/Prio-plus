#include "server.h"

#include <netinet/in.h>
#include <sys/socket.h>
#include <unistd.h>

#include <cstdlib>
#include <iostream>
#include <unordered_map>
#include <string>
#include <vector>
#include <queue>

#include "correlated.h"
#include "correlated_validated.h"
#include "eval_heavy.h"
#include "hash.h"
#include "heavy.h"
#include "net_share.h"
#include "ot.h"
#include "types.h"
#include "utils.h"

#define SERVER0_IP "127.0.0.1"
#define SERVER1_IP "127.0.0.1"

// #define SERVER0_IP "52.87.230.64"
// #define SERVER1_IP "54.213.189.18"

// Fail if more than this fraction of clients provide invalid inputs
#define INVALID_THRESHOLD 0.5

MultEvalManager* mult_eval_manager;

OT_Wrapper* ot0;
OT_Wrapper* ot1;
NetIO* garbleIO;

// Store wrapper for share conversion.
CorrelatedStore* correlated_store;
#define CACHE_SIZE 8192
// #define CACHE_SIZE 65536
// #define CACHE_SIZE 262144
// #define CACHE_SIZE 2097152
// If set, does fast but insecure offline precompute.
const int lazy_precompute = 0;

// Note: Currently does it in a single batch.
// I.e. receive and store all, then process all.
// TODO: Set up batching. Either every N inputs (based on space consumed),
// or every S seconds (figure out timer)

size_t send_tag(const int sockfd, const std::string x) {
    int ret = send(sockfd, &x[0], TAG_LENGTH, 0);
    if (ret <= 0) error_exit("Failed to send");
    return ret;
}

std::string recv_tag(const int serverfd) {
    char tag_buf[TAG_LENGTH];
    recv_in(serverfd, &tag_buf[0], TAG_LENGTH);
    std::string tag(tag_buf, tag_buf + TAG_LENGTH);
    return tag;
}

void bind_and_listen(sockaddr_in& addr, int& sockfd, const int port, const int reuse = 1) {
    sockfd = socket(AF_INET, SOCK_STREAM, 0);

    if (sockfd < 0) error_exit("Socket creation failed");

    if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, &reuse, sizeof(reuse)))
        error_exit("Sockopt failed");
    if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEPORT, &reuse, sizeof(reuse)))
        error_exit("Sockopt failed");

    bzero((char *) &addr, sizeof(addr));
    addr.sin_family = AF_INET;
    addr.sin_addr.s_addr = INADDR_ANY;
    addr.sin_port = htons(port);

    if (bind(sockfd, (sockaddr*)&addr, sizeof(addr)) < 0) {
        std::cerr << "Failed to bind to port: " << port << std::endl;
        error_exit("Bind to port failed");
    }

    if (listen(sockfd, 2) < 0)
        error_exit("Listen failed");
}

// Asymmetric: 1 connects to 0, 0 listens to 1.
void server0_listen(int& sockfd, int& newsockfd, const int port, const int reuse = 0) {
    sockaddr_in addr;
    bind_and_listen(addr, sockfd, port, reuse);

    socklen_t addrlen = sizeof(addr);
    std::cout << "  Waiting to accept\n";
    newsockfd = accept(sockfd, (sockaddr*)&addr, &addrlen);
    if (newsockfd < 0) error_exit("Accept failure");
    std::cout << "  Accepted\n";
}

void server1_connect(int& sockfd, const int port, const int reuse = 0) {
    sockfd = socket(AF_INET, SOCK_STREAM, 0);
    if (sockfd < 0) error_exit("Socket creation failed");

    if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, &reuse, sizeof(reuse)))
        error_exit("Sockopt failed");
    if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEPORT, &reuse, sizeof(reuse)))
        error_exit("Sockopt failed");

    sockaddr_in addr;
    bzero((char *) &addr, sizeof(addr));
    addr.sin_family = AF_INET;
    addr.sin_port = htons(port);
    inet_pton(AF_INET, SERVER0_IP, &addr.sin_addr);

    std::cout << "  Trying to connect...\n";
    if (connect(sockfd, (sockaddr*)&addr, sizeof(addr)) < 0)
        error_exit("Can't connect to other server");
    std::cout << "  Connected\n";
}

int recv_unvalidated(const int clientfd, const std::string tag, const size_t n) {
    if (STORE_TYPE != validate_store)
        return 0;

    // Overkill? Can be more efficient in the future. See client for more.
    const size_t N = NextPowerOfTwo(n-1);

    int cli_bytes = 0;

    DaBit** const bits = new DaBit*[N];
    AltTriple** const trips = new AltTriple*[N];
    for (unsigned int i = 0; i < N; i++) {
        bits[i] = new DaBit();
        trips[i] = new AltTriple();
    }

    cli_bytes += recv_DaBit_batch(clientfd, bits, N);
    cli_bytes += recv_AltTriple_batch(clientfd, trips, N);

    ((ValidateCorrelatedStore*) correlated_store)->queue_Unvalidated(bits, trips, tag, N);

    // std::cout << "got " << N << " unvalidated in " << cli_bytes << " bytes" << std::endl;

    return cli_bytes;
}

void process_unvalidated(const std::string tag, const size_t n) {
    if (STORE_TYPE != validate_store)
        return;

    ((ValidateCorrelatedStore*) correlated_store)->process_Unvalidated(tag);
}

// Batch of N (snips + num_input wire/share) validations
// Due to the nature of the final swap, both servers get the same valid array
// TODO: round split
const bool* const validate_snips(
        const size_t N, const size_t num_inputs,
        const int serverfd, const int server_num,
        Circuit* const * const circuit,
        const ClientPacket* const * const packet,
        const fmpz_t* const shares_p) {
    auto start = clock_start();

    // Setup precompute
    // TODO: support multiple circuits at once
    const size_t NumRoots = NextPowerOfTwo(circuit[0]->NumMulGates());
    mult_eval_manager->check_eval_point(3 * N);
    MultCheckPreComp* const pre = mult_eval_manager->get_Precomp(NumRoots);
    // Setup checker and cor
    Checker** const checker = new Checker*[N];
    Cor** const cor = new Cor*[N];
    for (unsigned int i = 0; i < N; i++) {
        checker[i] = new Checker(circuit[i], server_num, packet[i], pre,
                                 &shares_p[i * num_inputs]);
        cor[i] = checker[i]->CorFn();
    }

    // Reveal Cor
    reveal_Cor_batch(serverfd, cor, N);

    fmpz_t* valid_share; new_fmpz_array(&valid_share, N);
    for (unsigned int i = 0; i < N; i++) {
        checker[i]->OutShare(valid_share[i], cor[i]);
        delete cor[i];
        delete checker[i];
    }
    delete[] cor;
    delete[] checker;

    // Reveal valid
    reveal_fmpz_batch(serverfd, valid_share, N);

    // Answer is if valids are zero
    bool* const ans = new bool[N];
    for (unsigned int i = 0; i < N; i++) {
        ans[i] = fmpz_is_zero(valid_share[i]);
    }
    clear_fmpz_array(valid_share, N);

    std::cout << "snip circuit time: " << sec_from(start) << std::endl;
    return ans;
}

// N inputs sized M x L
// Accum M, so across N and L
// Vs accumulate, which is N sized M, with valid
void sum_accum(
        const size_t N, const size_t M, const size_t L,
        const fmpz_t* const x, fmpz_t* const accum) {
    for (unsigned int i = 0; i < N * M * L; i++)
        fmpz_mod_add(accum[(i/L)%M], accum[(i/L)%M], x[i], mod_ctx);
}

// Note: since bits, uses specific bitsum_ot, rather than normal store stuff.
returnType bit_sum(const initMsg msg, const int clientfd, const int serverfd,
                   const int server_num, uint64_t& ans) {
    std::unordered_map<std::string, bool> share_map;
    auto start = clock_start();

    BitShare share;
    const unsigned int total_inputs = msg.num_of_inputs;

    int cli_bytes = 0;
    for (unsigned int i = 0; i < total_inputs; i++) {
        cli_bytes += recv_in(clientfd, &share, sizeof(BitShare));
        const std::string tag(share.tag, share.tag + TAG_LENGTH);
        if (share_map.find(tag) != share_map.end())
            continue;
        share_map[tag] = share.val;
    }

    std::cout << "Received " << total_inputs << " total shares" << std::endl;
    std::cout << "bytes from client: " << cli_bytes << std::endl;
    std::cout << "receive time: " << sec_from(start) << std::endl;
    start = clock_start();
    auto start2 = clock_start();

    // if (STORE_TYPE != ot_store)
    //     ((CorrelatedStore*) correlated_store)->check_DaBits(total_inputs);

    int sent_bytes = 0;

    if (server_num == 1) {
        const size_t num_inputs = share_map.size();
        sent_bytes += send_size(serverfd, num_inputs);
        bool* const shares = new bool[num_inputs];
        int i = 0;
        for (const auto& share : share_map) {
            sent_bytes += send_tag(serverfd, share.first);
            shares[i] = share.second;
            i++;
        }
        std::cout << "tag time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        const uint64_t b = bitsum_ot_receiver(ot0, shares, num_inputs);
        delete[] shares;

        send_uint64(serverfd, b);
        std::cout << "convert time: " << sec_from(start2) << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;
        std::cout << "sent server bytes: " << sent_bytes << std::endl;
        return RET_NO_ANS;
    } else {
        size_t num_inputs, num_valid = 0;
        recv_size(serverfd, num_inputs);
        bool* const shares = new bool[num_inputs];
        bool* const valid = new bool[num_inputs];

        for (unsigned int i = 0; i < num_inputs; i++) {
            const std::string tag = recv_tag(serverfd);

            bool is_valid = (share_map.find(tag) != share_map.end());
            valid[i] = is_valid;
            if (!is_valid)
                continue;
            num_valid++;
            shares[i] = share_map[tag];
        }
        std::cout << "tag time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        const uint64_t a = bitsum_ot_sender(ot0, shares, valid, num_inputs);
        delete[] shares;
        delete[] valid;

        uint64_t b;
        recv_uint64(serverfd, b);

        std::cout << "Final valid count: " << num_valid << " / " << total_inputs << std::endl;
        std::cout << "convert time: " << sec_from(start2) << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;
        if (num_valid < total_inputs * (1 - INVALID_THRESHOLD)) {
            std::cout << "Failing, This is less than the invalid threshold of " << INVALID_THRESHOLD << std::endl;
            return RET_INVALID;
        }

        ans = a + b;
        return RET_ANS;
    }
}

returnType int_sum(const initMsg msg, const int clientfd, const int serverfd,
                   const int server_num, uint64_t& ans) {
    std::unordered_map<std::string, uint64_t> share_map;
    auto start = clock_start();

    IntShare share;
    const size_t num_bits = msg.num_bits;
    const uint64_t max_val = 1ULL << num_bits;
    const unsigned int total_inputs = msg.num_of_inputs;

    int cli_bytes = 0;
    for (unsigned int i = 0; i < total_inputs; i++) {
        cli_bytes += recv_in(clientfd, &share, sizeof(IntShare));
        const std::string tag(share.tag, share.tag + TAG_LENGTH);

        if (share_map.find(tag) != share_map.end()
            or share.val >= max_val)
            continue;
        share_map[tag] = share.val;

        // std::cout << "share[" << i << "] = " << share.val << std::endl;

        cli_bytes += recv_unvalidated(clientfd, tag, num_bits);
    }

    std::cout << "Received " << total_inputs << " total shares" << std::endl;
    std::cout << "bytes from client: " << cli_bytes << std::endl;
    std::cout << "receive time: " << sec_from(start) << std::endl;

    if (STORE_TYPE == precompute_store)
        ((CorrelatedStore*) correlated_store)->check_DaBits(total_inputs * num_bits);

    std::cout << "Sizes after recv:" << std::endl;
    correlated_store->print_Sizes();

    start = clock_start();
    auto start2 = clock_start();
    int sent_bytes = 0;
    size_t num_inputs;

    if (server_num == 1) {
        num_inputs = share_map.size();
        sent_bytes += send_size(serverfd, num_inputs);
    } else {
        recv_size(serverfd, num_inputs);
    }
    uint64_t* const shares = new uint64_t[num_inputs];
    bool* const valid = new bool[num_inputs];
    memset(valid, true, num_inputs * sizeof(bool));

    /* Share Sync */
    if (server_num == 1) {
        int i = 0;
        for (const auto& share : share_map) {
            sent_bytes += send_tag(serverfd, share.first);
            shares[i] = share.second;
            i++;

            process_unvalidated(&share.first[0], num_bits);
        }
    } else {
        for (unsigned int i = 0; i < num_inputs; i++) {
            const std::string tag = recv_tag(serverfd);

            bool is_valid = (share_map.find(tag) != share_map.end());
            valid[i] = is_valid;
            if (!is_valid)
                continue;
            shares[i] = share_map[tag];

            process_unvalidated(tag, num_bits);
        }
    }
    std::cout << "tag time: " << sec_from(start2) << std::endl;
    start2 = clock_start();

    fmpz_t* shares_p; new_fmpz_array(&shares_p, num_inputs);
    fmpz_t* flat_xp; new_fmpz_array(&flat_xp, num_inputs * num_bits);
    bool* const send_buff = new bool[2 * num_inputs * num_bits];
    bool* const recv_buff = new bool[2 * num_inputs * num_bits];
    memset(recv_buff, 0, 2 * num_inputs * num_bits);

    if (STORE_TYPE == validate_store) {
        // TODO: batch validation in parallel with conversion
        ValidateCorrelatedStore* s = (ValidateCorrelatedStore*) correlated_store;

        s->batch_Validate(num_inputs * num_bits);

        correlated_store->b2a_multi_setup(num_inputs, num_bits, shares, flat_xp, send_buff);
        memcpy(&send_buff[num_inputs], valid, num_inputs);

        sent_bytes += swap_bool_batch(serverfd, send_buff, recv_buff, 2 * num_inputs);

        correlated_store->b2a_multi_finish(num_inputs, num_bits, shares_p, flat_xp, send_buff, recv_buff);
        delete[] send_buff;
        clear_fmpz_array(flat_xp, num_inputs * num_bits);
        for (unsigned int i = 0; i < num_inputs; i++)
            valid[i] &= recv_buff[num_inputs + i];
        delete[] recv_buff;
    } else {
        correlated_store->b2a_multi_setup(num_inputs, num_bits, shares, flat_xp, send_buff);
        memcpy(&send_buff[num_inputs], valid, num_inputs);

        sent_bytes += swap_bool_batch(serverfd, send_buff, recv_buff, 2 * num_inputs);

        correlated_store->b2a_multi_finish(num_inputs, num_bits, shares_p, flat_xp, send_buff, recv_buff);
        delete[] send_buff;
        clear_fmpz_array(flat_xp, num_inputs * num_bits);
        for (unsigned int i = 0; i < num_inputs; i++)
            valid[i] &= recv_buff[num_inputs + i];
        delete[] recv_buff;
    }

    // Then marge in batch validation (two rounds)
    sent_bytes += correlated_store->b2a_multi(num_inputs, num_bits, shares, shares_p);

    std::cout << "convert time: " << sec_from(start2) << std::endl;
    start2 = clock_start();

    fmpz_t* b; new_fmpz_array(&b, 1);
    int num_valid = accumulate(num_inputs, 1, shares_p, valid, b);
    clear_fmpz_array(shares_p, num_inputs);
    delete[] shares;
    delete[] valid;
    std::cout << "accumulate time: " << sec_from(start2) << std::endl;

    std::cout << "Final valid count: " << num_valid << " / " << total_inputs << std::endl;
    std::cout << "total compute time: " << sec_from(start) << std::endl;
    std::cout << "sent server bytes: " << sent_bytes << std::endl;
    if (num_valid < total_inputs * (1 - INVALID_THRESHOLD)) {
        std::cout << "Failing, This is less than the invalid threshold of " << INVALID_THRESHOLD << std::endl;
        clear_fmpz_array(b, 1);
        return RET_INVALID;
    }

    reveal_fmpz_batch(serverfd, b, 1);
    ans = fmpz_get_ui(b[0]);
    clear_fmpz_array(b, 1);
    return RET_ANS;
}

// For AND and OR
returnType xor_op(const initMsg msg, const int clientfd, const int serverfd,
                  const int server_num, bool& ans) {
    std::unordered_map<std::string, uint64_t> share_map;
    auto start = clock_start();

    IntShare share;
    const unsigned int total_inputs = msg.num_of_inputs;

    int cli_bytes = 0;
    for (unsigned int i = 0; i < total_inputs; i++) {
        cli_bytes += recv_in(clientfd, &share, sizeof(IntShare));
        const std::string tag(share.tag, share.tag + TAG_LENGTH);

        if (share_map.find(tag) != share_map.end())
            continue;
        share_map[tag] = share.val;
    }

    std::cout << "Received " << total_inputs << " total shares" << std::endl;
    std::cout << "bytes from client: " << cli_bytes << std::endl;
    std::cout << "receive time: " << sec_from(start) << std::endl;
    start = clock_start();
    auto start2 = clock_start();

    int sent_bytes = 0;

    if (server_num == 1) {
        const size_t num_inputs = share_map.size();
        sent_bytes += send_size(serverfd, num_inputs);
        uint64_t b = 0;
        uint64_t* const shares = new uint64_t[num_inputs];
        size_t idx = 0;
        for (const auto& share : share_map) {
            sent_bytes += send_tag(serverfd, share.first);
            shares[idx] = share.second;
            idx++;
        }
        std::cout << "tag time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        bool* const other_valid = new bool[num_inputs];
        recv_bool_batch(serverfd, other_valid, num_inputs);
        for (unsigned int i = 0; i < num_inputs; i++) {
            if (!other_valid[i])
                continue;
            b ^= shares[i];
        }
        delete[] other_valid;
        delete[] shares;

        send_uint64(serverfd, b);
        std::cout << "convert time: " << sec_from(start2) << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;
        std::cout << "sent server bytes: " << sent_bytes << std::endl;
        return RET_NO_ANS;
    } else {
        size_t num_inputs, num_valid = 0;
        recv_size(serverfd, num_inputs);
        uint64_t a = 0;
        bool* const valid = new bool[num_inputs];

        for (unsigned int i = 0; i < num_inputs; i++) {
            const std::string tag = recv_tag(serverfd);
            valid[i] = (share_map.find(tag) != share_map.end());
            if (!valid[i])
                continue;
            num_valid++;
            a ^= share_map[tag];
        }

        std::cout << "tag + convert time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        sent_bytes += send_bool_batch(serverfd, valid, num_inputs);

        delete[] valid;

        uint64_t b;
        recv_uint64(serverfd, b);
        const uint64_t aggr = a ^ b;

        std::cout << "Final valid count: " << num_valid << " / " << total_inputs << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;
        std::cout << "sent server bytes: " << sent_bytes << std::endl;
        if (num_valid < total_inputs * (1 - INVALID_THRESHOLD)) {
            std::cout << "Failing, This is less than the invalid threshold of " << INVALID_THRESHOLD << std::endl;
            return RET_INVALID;
        }

        if (msg.type == AND_OP) {
            ans = (aggr == 0);
        } else if (msg.type == OR_OP) {
            ans = (aggr != 0);
        } else {
            error_exit("Message type incorrect for xor_op");
        }
        return RET_ANS;
    }
}

// For MAX and MIN
returnType max_op(const initMsg msg, const int clientfd, const int serverfd,
                  const int server_num, uint64_t& ans) {
    std::unordered_map<std::string, uint64_t*> share_map;
    auto start = clock_start();

    MaxShare share;
    const unsigned int total_inputs = msg.num_of_inputs;
    const unsigned int B = msg.max_inp;
    // Need this to have all share arrays stay in memory, for server1 later.
    uint64_t* const shares = new uint64_t[total_inputs * (B + 1)];

    int cli_bytes = 0;
    for (unsigned int i = 0; i < total_inputs; i++) {
        cli_bytes += recv_in(clientfd, &share.tag[0], TAG_LENGTH);
        const std::string tag(share.tag, share.tag + TAG_LENGTH);

        cli_bytes += recv_uint64_batch(clientfd, &shares[i*(B+1)], B+1);

        if (share_map.find(tag) != share_map.end())
            continue;
        share_map[tag] = &shares[i*(B+1)];
    }

    std::cout << "Received " << total_inputs << " total shares" << std::endl;
    std::cout << "bytes from client: " << cli_bytes << std::endl;
    std::cout << "receive time: " << sec_from(start) << std::endl;
    start = clock_start();
    auto start2 = clock_start();

    int sent_bytes = 0;

    if (server_num == 1) {
        const size_t num_inputs = share_map.size();
        sent_bytes += send_size(serverfd, num_inputs);
        uint64_t b[B+1];
        memset(b, 0, sizeof(b));
        uint64_t** const ordered_shares = new uint64_t*[num_inputs];
        size_t idx = 0;
        for (const auto& share : share_map) {
            sent_bytes += send_tag(serverfd, share.first);
            ordered_shares[idx] = share.second;
            idx++;
        }
        std::cout << "tag time: " << sec_from(start2) << std::endl;
        start2 = clock_start();
        bool* const other_valid = new bool[num_inputs];
        recv_bool_batch(serverfd, other_valid, num_inputs);
        for (unsigned int i = 0; i < num_inputs; i++) {
            if (!other_valid[i])
                continue;
            for (unsigned int j = 0; j <= B; j++)
                b[j] ^= ordered_shares[i][j];
        }
        delete[] other_valid;
        delete[] shares;
        delete[] ordered_shares;
        send_uint64_batch(serverfd, b, B+1);

        std::cout << "convert time: " << sec_from(start2) << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;
        std::cout << "sent server bytes: " << sent_bytes << std::endl;
        return RET_NO_ANS;
    } else {
        size_t num_inputs, num_valid = 0;
        recv_size(serverfd, num_inputs);
        uint64_t a[B+1];
        memset(a, 0, sizeof(a));
        bool* const valid = new bool[num_inputs];

        for (unsigned int i = 0; i < num_inputs; i++) {
            const std::string tag = recv_tag(serverfd);
            valid[i] = (share_map.find(tag) != share_map.end());
            if (!valid[i])
                continue;
            num_valid++;
            for (unsigned int j = 0; j <= B; j++)
                a[j] ^= share_map[tag][j];
        }

        std::cout << "tag+convert time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        sent_bytes += send_bool_batch(serverfd, valid, num_inputs);

        delete[] shares;
        delete[] valid;
        uint64_t b[B+1];
        recv_uint64_batch(serverfd, b, B+1);

        std::cout << "Final valid count: " << num_valid << " / " << total_inputs << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;
        std::cout << "sent server bytes: " << sent_bytes << std::endl;
        if (num_valid < total_inputs * (1 - INVALID_THRESHOLD)) {
            std::cout << "Failing, This is less than the invalid threshold of " << INVALID_THRESHOLD << std::endl;
            return RET_INVALID;
        }

        for (unsigned int j = B; j >= 0; j--) {
            // std::cout << "a, b[" << j << "] = " << a[j] << ", " << b[j] << std::endl;
            if (a[j] != b[j]) {
                if (msg.type == MAX_OP) {
                    ans = j;
                } else if (msg.type == MIN_OP) {
                    ans = B - j;
                } else {
                    error_exit("Message type incorrect for max_op");
                }
                return RET_ANS;
            }
        }
        // Shouldn't reach.
        return RET_INVALID;
    }

    return RET_INVALID;
}

// For var, stddev
returnType var_op(const initMsg msg, const int clientfd, const int serverfd,
                  const int server_num, double& ans) {
    auto start = clock_start();

    typedef std::tuple <uint64_t, uint64_t, ClientPacket*> sharetype;
    std::unordered_map<std::string, sharetype> share_map;

    VarShare share;
    const size_t num_bits = msg.num_bits;
    const uint64_t max_val = 1ULL << num_bits;
    const unsigned int total_inputs = msg.num_of_inputs;
    const size_t num_dabits = 3;

    // Just for getting sizes
    const Circuit* const mock_circuit = CheckVar();
    const size_t NMul = mock_circuit->NumMulGates();
    delete mock_circuit;

    int cli_bytes = 0;
    for (unsigned int i = 0; i < total_inputs; i++) {
        cli_bytes += recv_in(clientfd, &share, sizeof(VarShare));
        const std::string tag(share.tag, share.tag + TAG_LENGTH);

        ClientPacket* const packet = new ClientPacket(NMul);
        int packet_bytes = recv_ClientPacket(clientfd, packet, NMul);
        cli_bytes += packet_bytes;

        if ((share_map.find(tag) != share_map.end())
            or (share.val >= max_val)
            or (share.val_squared >= max_val * max_val)
            or (packet_bytes <= 0)
            ) {
            delete packet;
            continue;
        }
        share_map[tag] = {share.val, share.val_squared, packet};

        cli_bytes += recv_unvalidated(clientfd, tag, num_bits * num_dabits);
    }

    std::cout << "Received " << total_inputs << " total shares" << std::endl;
    std::cout << "bytes from client: " << cli_bytes << std::endl;
    std::cout << "receive time: " << sec_from(start) << std::endl;

    if (STORE_TYPE != ot_store)
        ((CorrelatedStore*) correlated_store)->check_DaBits(total_inputs * num_bits * num_dabits);

    start = clock_start();
    auto start2 = clock_start();
    int sent_bytes = 0;
    size_t num_inputs;

    if (server_num == 1) {
        num_inputs = share_map.size();
        sent_bytes += send_size(serverfd, num_inputs);
    } else {
        recv_size(serverfd, num_inputs);
    }

    uint64_t* const shares = new uint64_t[2 * num_inputs];
    bool* const valid = new bool[num_inputs];
    memset(valid, true, num_inputs * sizeof(bool));
    ClientPacket** const packet = new ClientPacket*[num_inputs];
    Circuit** const circuit = new Circuit*[num_inputs];

    if (server_num == 1) {
        size_t i = 0;
        for (const auto& share : share_map) {
            sent_bytes += send_tag(serverfd, share.first);

            std::tie(shares[i], shares[num_inputs + i], packet[i]) = share.second;
            circuit[i] = CheckVar();

            i++;

            process_unvalidated(&share.first[0], num_bits * num_dabits);
        }
    } else {
        for (unsigned int i = 0; i < num_inputs; i++) {
            const std::string tag = recv_tag(serverfd);

            valid[i] = (share_map.find(tag) != share_map.end());
            if (valid[i]) {
                std::tie(shares[i], shares[num_inputs + i], packet[i]) = share_map[tag];
            } else {
                packet[i] = new ClientPacket(NMul);  // mock empty packet
            }
            circuit[i] = CheckVar();

            process_unvalidated(tag, num_bits * num_dabits);
        }
    }
    std::cout << "tag time: " << sec_from(start2) << std::endl;
    start2 = clock_start();

    fmpz_t* flat_xp_1; new_fmpz_array(&flat_xp_1, num_inputs * num_bits);
    bool* const v = new bool[num_inputs * num_bits * 3];
    correlated_store->b2a_multi_setup(num_inputs, num_bits, shares, flat_xp_1, v);
    fmpz_t* flat_xp_2; new_fmpz_array(&flat_xp_2, num_inputs * num_bits * 2);
    correlated_store->b2a_multi_setup(num_inputs, 2 * num_bits, &shares[num_inputs],
            flat_xp_2, &v[num_inputs * num_bits]);
    delete[] shares;
    bool* const v_other = new bool[num_inputs * num_bits * 3];

    sent_bytes += swap_bool_batch(serverfd, v, v_other, num_inputs * num_bits * 3);

    fmpz_t* shares_p; new_fmpz_array(&shares_p, num_inputs * 2);
    correlated_store->b2a_multi_finish(num_inputs, num_bits, shares_p, flat_xp_1, v, v_other);
    clear_fmpz_array(flat_xp_1, num_inputs * num_bits);
    correlated_store->b2a_multi_finish(num_inputs, 2 * num_bits, &shares_p[num_inputs],
            flat_xp_2, &v[num_inputs * num_bits], &v_other[num_inputs * num_bits]);
    clear_fmpz_array(flat_xp_2, num_inputs * num_bits * 2);
    delete[] v;
    delete[] v_other;
    transpose_fmpz_array(shares_p, num_inputs, 2);
    std::cout << "convert time: " << sec_from(start2) << std::endl;
    start2 = clock_start();

    const bool* const snip_valid = validate_snips(
        num_inputs, 2, serverfd, server_num, circuit, packet, shares_p);

    for (unsigned int i = 0; i < num_inputs; i++) {
        delete circuit[i];
        delete packet[i];
    }
    delete[] circuit;
    delete[] packet;

    if (server_num == 1) {
        recv_bool_batch(serverfd, valid, num_inputs);

        for (unsigned int i = 0; i < num_inputs; i++)
            valid[i] &= snip_valid[i];
    } else {
        for (unsigned int i = 0; i < num_inputs; i++)
            valid[i] &= snip_valid[i];

        sent_bytes += send_bool_batch(serverfd, valid, num_inputs);
    }
    delete[] snip_valid;
    std::cout << "validate time: " << sec_from(start2) << std::endl;
    start2 = clock_start();

    fmpz_t* b; new_fmpz_array(&b, 2);
    const size_t num_valid = accumulate(num_inputs, 2, shares_p, valid, b);
    delete[] valid;
    clear_fmpz_array(shares_p, num_inputs * 2);
    std::cout << "accumulate time: " << sec_from(start2) << std::endl;
    start2 = clock_start();

    std::cout << "Final valid count: " << num_valid << " / " << total_inputs << std::endl;
    std::cout << "total compute time: " << sec_from(start) << std::endl;
    std::cout << "sent non-snip server bytes: " << sent_bytes << std::endl;
    if (num_valid < total_inputs * (1 - INVALID_THRESHOLD)) {
        std::cout << "Failing, This is less than the invalid threshold of " << INVALID_THRESHOLD << std::endl;
        clear_fmpz_array(b, 2);
        return RET_INVALID;
    }

    reveal_fmpz_batch(serverfd, b, 2);

    const double ex = fmpz_get_d(b[0]) / num_valid;
    const double ex2 = fmpz_get_d(b[1]) / num_valid;
    ans = ex2 - (ex * ex);
    if (msg.type == VAR_OP) {
        std::cout << "Ans: " << ex2 << " - (" << ex << ")^2 = " << ans << std::endl;
    }
    if (msg.type == STDDEV_OP) {
        ans = sqrt(ans);
        std::cout << "Ans: sqrt(" << ex2 << " - (" << ex << ")^2) = " << ans << std::endl;
    }
    clear_fmpz_array(b, 2);
    std::cout << "eval time: " << sec_from(start2) << std::endl;
    return RET_ANS;
}

returnType linreg_op(const initMsg msg, const int clientfd,
                     const int serverfd, const int server_num) {
    auto start = clock_start();
    int cli_bytes = 0;

    size_t degree;
    cli_bytes += recv_size(clientfd, degree);

    std::cout << "Linreg degree: " << degree << std::endl;

    const size_t num_x = degree - 1;
    const size_t num_quad = num_x * (num_x + 1) / 2;
    const size_t num_fields = 2 * num_x + 1 + num_quad;
    // std::cout << "num_x: " << num_x << std::endl;
    // std::cout << "num_quad: " << num_quad << std::endl;
    // std::cout << "num_fields: " << num_fields << std::endl;
    const int num_dabits = degree * (degree + 2) - 2;

    // [x], y, [x2], [xy]
    typedef std::tuple <uint64_t*, uint64_t, uint64_t*, uint64_t*, ClientPacket*> sharetype;
    std::unordered_map<std::string, sharetype> share_map;

    const size_t num_bits = msg.num_bits;
    const uint64_t max_val = 1ULL << num_bits;
    const unsigned int total_inputs = msg.num_of_inputs;

    // Just for getting sizes
    Circuit* const mock_circuit = CheckLinReg(degree);
    const size_t NMul = mock_circuit->NumMulGates();
    delete mock_circuit;

    LinRegShare share;
    for (unsigned int i = 0; i < total_inputs; i++) {
        bool sizes_valid = true;

        cli_bytes += recv_in(clientfd, &share.tag[0], TAG_LENGTH);
        const std::string tag(share.tag, share.tag + TAG_LENGTH);

        share.x_vals = new uint64_t[num_x];
        share.x2_vals = new uint64_t[num_quad];
        share.xy_vals = new uint64_t[num_x];

        cli_bytes += recv_uint64_batch(clientfd, share.x_vals, num_x);
        cli_bytes += recv_uint64(clientfd, share.y);
        cli_bytes += recv_uint64_batch(clientfd, share.x2_vals, num_quad);
        cli_bytes += recv_uint64_batch(clientfd, share.xy_vals, num_x);

        for (unsigned int j = 0; j < num_x; j++) {
            if (share.x_vals[j] >= max_val)
                sizes_valid = false;
            if (share.xy_vals[j] >= max_val * max_val)
                sizes_valid = false;
        }
        if (share.y >= max_val)
            sizes_valid = false;
        for (unsigned int j = 0; j < num_quad; j++) {
            if (share.x2_vals[j] >= max_val * max_val)
                sizes_valid = false;
        }

        ClientPacket* const packet = new ClientPacket(NMul);
        int packet_bytes = recv_ClientPacket(clientfd, packet, NMul);
        cli_bytes += packet_bytes;

        if ((share_map.find(tag) != share_map.end())
            or (not sizes_valid)
            or (packet_bytes  <= 0)
            ) {
            delete[] share.x_vals;
            delete[] share.x2_vals;
            delete[] share.xy_vals;
            delete packet;
            continue;
        }

        share_map[tag] = {share.x_vals, share.y, share.x2_vals, share.xy_vals, packet};

        cli_bytes += recv_unvalidated(clientfd, tag, num_bits * num_dabits);
    }

    std::cout << "Received " << total_inputs << " total shares" << std::endl;
    std::cout << "bytes from client: " << cli_bytes << std::endl;
    std::cout << "receive time: " << sec_from(start) << std::endl;

    if (STORE_TYPE != ot_store)
        ((CorrelatedStore*) correlated_store)->check_DaBits(total_inputs * num_bits * num_dabits);

    start = clock_start();
    auto start2 = clock_start();
    int sent_bytes = 0;
    size_t num_inputs;

    if (server_num == 1) {
        num_inputs = share_map.size();
        sent_bytes += send_size(serverfd, num_inputs);
    } else {
        recv_size(serverfd, num_inputs);
    }

    uint64_t* const shares = new uint64_t[num_inputs * num_fields];
    bool* const valid = new bool[num_inputs];
    memset(valid, true, num_inputs * sizeof(bool));
    ClientPacket** const packet = new ClientPacket*[num_inputs];
    Circuit** const circuit = new Circuit*[num_inputs];

    if (server_num == 1) {
        size_t i = 0;
        for (const auto& share : share_map) {
            sent_bytes += send_tag(serverfd, share.first);

            circuit[i] = CheckLinReg(degree);

            memcpy(&shares[num_fields * i],
                   std::get<0>(share.second), num_x * sizeof(uint64_t));
            shares[num_fields * i + num_x] = std::get<1>(share.second);
            memcpy(&shares[num_fields * i + num_x + 1],
                   std::get<2>(share.second), num_quad * sizeof(uint64_t));
            memcpy(&shares[num_fields * i + num_x + num_quad + 1],
                   std::get<3>(share.second), num_x * sizeof(uint64_t));
            packet[i] = std::get<4>(share.second);

            i++;

            delete std::get<0>(share.second);
            delete std::get<2>(share.second);
            delete std::get<3>(share.second);

            process_unvalidated(&share.first[0], num_bits * num_dabits);
        }
    } else {
        for (unsigned int i = 0; i < num_inputs; i++) {
            const std::string tag = recv_tag(serverfd);

            valid[i] = (share_map.find(tag) != share_map.end());
            if (valid[i]) {
                memcpy(&shares[num_fields * i],
                   std::get<0>(share_map[tag]), num_x * sizeof(uint64_t));
                shares[num_fields * i + num_x] = std::get<1>(share_map[tag]);
                memcpy(&shares[num_fields * i + num_x + 1],
                       std::get<2>(share_map[tag]), num_quad * sizeof(uint64_t));
                memcpy(&shares[num_fields * i + num_x + num_quad + 1],
                       std::get<3>(share_map[tag]), num_x * sizeof(uint64_t));
                packet[i] = std::get<4>(share_map[tag]);

                delete std::get<0>(share_map[tag]);
                delete std::get<2>(share_map[tag]);
                delete std::get<3>(share_map[tag]);
            } else {
                packet[i] = new ClientPacket(NMul);
            }
            circuit[i] = CheckLinReg(degree);

            process_unvalidated(tag, num_bits * num_dabits);
        }
    }

    std::cout << "tag time: " << sec_from(start2) << std::endl;
    start2 = clock_start();

    const size_t num_1 = degree;
    const size_t num_2 = num_fields - num_1;
    const size_t total = num_bits * (num_1 + 2 * num_2);
    uint64_t* const shares_transposed = new uint64_t[num_inputs * num_fields];
    for (unsigned int i = 0; i < num_inputs; i++) {
        for (unsigned int j = 0; j < num_fields; j++) {
            shares_transposed[j * num_inputs + i] = shares[i * num_fields + j];
        }
    }
    delete[] shares;
    fmpz_t* flat_xp_1; new_fmpz_array(&flat_xp_1, num_inputs * num_1 * num_bits);
    bool* const v = new bool[num_inputs * total];
    correlated_store->b2a_multi_setup(num_inputs * num_1, num_bits, shares_transposed, flat_xp_1, v);
    fmpz_t* flat_xp_2; new_fmpz_array(&flat_xp_2, num_inputs * num_2 * 2 * num_bits);
    correlated_store->b2a_multi_setup(num_inputs * num_2, 2 * num_bits,
        &shares_transposed[num_inputs * num_1], flat_xp_2, &v[num_inputs * num_1 * num_bits]);
    // delete[] shares_transposed;
    bool* const v_other = new bool[num_inputs * total];

    sent_bytes += swap_bool_batch(serverfd, v, v_other, num_inputs * total);

    fmpz_t* shares_p; new_fmpz_array(&shares_p, num_inputs * num_fields);
    correlated_store->b2a_multi_finish(num_inputs * num_1, num_bits, shares_p, flat_xp_1, v, v_other);
    clear_fmpz_array(flat_xp_1, num_inputs * num_1 * num_bits);
    correlated_store->b2a_multi_finish(num_inputs * num_2, 2 * num_bits, &shares_p[num_inputs * num_1],
            flat_xp_2, &v[num_inputs * num_1 * num_bits], &v_other[num_inputs * num_1 * num_bits]);
    clear_fmpz_array(flat_xp_2, num_inputs * num_2 * 2 * num_bits);
    delete[] v;
    delete[] v_other;
    transpose_fmpz_array(shares_p, num_inputs, num_fields);

    std::cout << "convert time: " << sec_from(start2) << std::endl;
    start2 = clock_start();

    const bool* const snip_valid = validate_snips(
            num_inputs, num_fields, serverfd, server_num, circuit,
            packet, shares_p);

    for (unsigned int i = 0; i < num_inputs; i++) {
        delete circuit[i];
        delete packet[i];
    }
    delete[] circuit;
    delete[] packet;

    if (server_num == 1) {
        recv_bool_batch(serverfd, valid, num_inputs);

        for (unsigned int i = 0; i < num_inputs; i++)
            valid[i] &= snip_valid[i];
    } else {
        for (unsigned int i = 0; i < num_inputs; i++)
            valid[i] &= snip_valid[i];

        sent_bytes += send_bool_batch(serverfd, valid, num_inputs);
    }
    delete[] snip_valid;
    std::cout << "validate time: " << sec_from(start2) << std::endl;
    start2 = clock_start();

    fmpz_t* b; new_fmpz_array(&b, num_fields);
    size_t num_valid = accumulate(num_inputs, num_fields, shares_p, valid, b);
    delete[] valid;
    clear_fmpz_array(shares_p, num_inputs * num_fields);
    std::cout << "accumulate time: " << sec_from(start2) << std::endl;
    start2 = clock_start();

    std::cout << "Final valid count: " << num_valid << " / " << total_inputs << std::endl;
    std::cout << "total compute time: " << sec_from(start) << std::endl;
    std::cout << "sent non-snip server bytes: " << sent_bytes << std::endl;
    if (num_valid < total_inputs * (1 - INVALID_THRESHOLD)) {
        std::cout << "Failing, This is less than the invalid threshold of " << INVALID_THRESHOLD << std::endl;
        clear_fmpz_array(b, num_fields);
        return RET_INVALID;
    }

    reveal_fmpz_batch(serverfd, b, num_fields);

    uint64_t* const x_accum = new uint64_t[degree + num_quad];
    memset(x_accum, 0, (degree + num_quad) * sizeof(uint64_t));
    uint64_t* const y_accum = new uint64_t[degree];
    memset(y_accum, 0, (degree) * sizeof(uint64_t));

    x_accum[0] = num_valid;
    for (unsigned int j = 0; j < num_x; j++)
        x_accum[1 + j] = fmpz_get_si(b[j]);
    y_accum[0] = fmpz_get_si(b[num_x]);
    for (unsigned int j = 0; j < num_quad; j++)
        x_accum[1 + num_x + j] = fmpz_get_si(b[degree + j]);
    for (unsigned int j = 0; j < num_x; j++)
        y_accum[1 + j] = fmpz_get_si(b[degree + num_quad + j]);
    clear_fmpz_array(b, num_fields);

    const double* const c = SolveLinReg(degree, x_accum, y_accum);
    std::cout << "Estimate: y = ";
    for (unsigned int i = 0; i < degree; i++) {
        if (i > 0) std::cout << " + ";
        std::cout << c[i];
        if (i > 0) std::cout << " * x_" << (i-1);
    }
    std::cout << std::endl;
    std::cout << "eval time: " << sec_from(start2) << std::endl;
    delete[] x_accum;
    delete[] y_accum;
    delete[] c;

    return RET_ANS;
}

returnType freq_op(const initMsg msg, const int clientfd, const int serverfd,
                   const int server_num) {
    std::unordered_map<std::string, bool*> share_map;
    auto start = clock_start();

    const unsigned int total_inputs = msg.num_of_inputs;
    // const uint64_t max_inp = msg.max_inp;
    const uint64_t max_inp = 1ULL << msg.num_bits;

    FreqShare share;
    int cli_bytes = 0;
    for (unsigned int i = 0; i < total_inputs; i++) {
        cli_bytes += recv_in(clientfd, &share.tag[0], TAG_LENGTH);
        const std::string tag(share.tag, share.tag+TAG_LENGTH);
        share.arr = new bool[max_inp];
        cli_bytes += recv_bool_batch(clientfd, share.arr, max_inp);

        if (share_map.find(tag) != share_map.end()) {
            delete[] share.arr;
            continue;
        }
        share_map[tag] = share.arr;

        // for (unsigned int j = 0; j < max_inp; j++) {
        //     std::cout << "share[" << i << ", " << j << "] = " << share.arr[j] << std::endl;
        // }
    }

    std::cout << "Received " << total_inputs << " total shares" << std::endl;
    std::cout << "bytes from client: " << cli_bytes << std::endl;
    std::cout << "receive time: " << sec_from(start) << std::endl;

    if (STORE_TYPE != ot_store)
        ((CorrelatedStore*) correlated_store)->check_DaBits(total_inputs * max_inp);

    start = clock_start();
    auto start2 = clock_start();
    int sent_bytes = 0;

    if (server_num == 1) {
        const size_t num_inputs = share_map.size();
        sent_bytes += send_size(serverfd, num_inputs);
        bool* const shares = new bool[num_inputs * max_inp];

        size_t idx = 0;
        for (const auto& share : share_map) {
            sent_bytes += send_tag(serverfd, share.first);
            memcpy(&shares[idx * max_inp], share.second, max_inp);
            delete[] share.second;
            idx++;
        }
        std::cout << "tag time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        fmpz_t* shares_p; new_fmpz_array(&shares_p, num_inputs * max_inp);
        sent_bytes += correlated_store->b2a_single(num_inputs * max_inp, shares, shares_p);
        std::cout << "convert time: " << sec_from(start2) << std::endl;
        delete[] shares;
        start2 = clock_start();

        bool* const valid = new bool[num_inputs];
        fmpz_t* sums; new_fmpz_array(&sums, num_inputs);
        fmpz_t sum; fmpz_init_set_si(sum, 0);
        for (unsigned int i = 0; i < num_inputs; i++) {
            fmpz_zero(sums[i]);
            for (unsigned int j = 0; j < max_inp; j++) {
                fmpz_mod_add(sums[i], sums[i], shares_p[i * max_inp + j], mod_ctx);
            }
            fmpz_mod_add(sum, sum, sums[i], mod_ctx);
        }
        // Batch check.
        fmpz_t sum_other; fmpz_init(sum_other);
        sent_bytes += send_fmpz(serverfd, sum);
        recv_fmpz(serverfd, sum_other);

        fmpz_mod_add(sum, sum, sum_other, mod_ctx);
        bool all_valid = fmpz_equal_ui(sum, total_inputs);
        fmpz_clear(sum);
        if (!all_valid) {
            sent_bytes += send_fmpz_batch(serverfd, sums, num_inputs);
        }
        recv_bool_batch(serverfd, valid, num_inputs);

        clear_fmpz_array(sums, num_inputs);
        std::cout << "validate time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        fmpz_t* b; new_fmpz_array(&b, max_inp);
        accumulate(num_inputs, max_inp, shares_p, valid, b);
        std::cout << "accumulate time: " << sec_from(start2) << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;

        delete[] valid;
        clear_fmpz_array(shares_p, num_inputs * max_inp);
        // send b
        send_fmpz_batch(serverfd, b, max_inp);
        clear_fmpz_array(b, max_inp);

        std::cout << "sent server bytes: " << sent_bytes << std::endl;
        return RET_NO_ANS;
    } else {
        size_t num_inputs;
        recv_size(serverfd, num_inputs);
        bool* const shares = new bool[num_inputs * max_inp];
        bool* const valid = new bool[num_inputs];

        for (unsigned int i = 0; i < num_inputs; i++) {
            const std::string tag = recv_tag(serverfd);
            valid[i] = (share_map.find(tag) != share_map.end());

            // realign shares_2 to tag order
            if (valid[i]) {
                memcpy(&shares[i * max_inp], share_map[tag], max_inp);
                delete[] share_map[tag];
            } else {
                memset(&shares[i * max_inp], 0, max_inp);
            }
        }
        std::cout << "tag time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        fmpz_t* shares_p; new_fmpz_array(&shares_p, num_inputs * max_inp);
        sent_bytes += correlated_store->b2a_single(num_inputs * max_inp, shares, shares_p);
        std::cout << "convert time: " << sec_from(start2) << std::endl;
        delete[] shares;
        start2 = clock_start();


        /* N freq vectors, each length M
           Check that arithmetic shares  sum to 1

           Originally:
            - Also check that boolean shares have odd parity
            - not actually be needed, as long as M < Int_modulus
            - just stops one wraparound.

           Here: two layer check
             if usually valid, check 1 is global agg.
             Check 2 is individual
             check 2 can also be binary search, but a lot more rounds, and bad worst case
             Becomes harder to integrate into round overlap (conditional rounds)
             worse if any bad expected
             Requires N * M < Int_modulus, to avoid rollover
        */
        // Validity check
        fmpz_t* sums; new_fmpz_array(&sums, num_inputs);
        fmpz_t sum; fmpz_init_set_si(sum, 0);
        for (unsigned int i = 0; i < num_inputs; i++) {
            fmpz_zero(sums[i]);
            for (unsigned int j = 0; j < max_inp; j++) {
                fmpz_mod_add(sums[i], sums[i], shares_p[i * max_inp + j], mod_ctx);
            }
            fmpz_mod_add(sum, sum, sums[i], mod_ctx);
        }
        // batch check
        fmpz_t sum_other; fmpz_init(sum_other);
        sent_bytes += send_fmpz(serverfd, sum);
        recv_fmpz(serverfd, sum_other);

        fmpz_mod_add(sum, sum, sum_other, mod_ctx);
        bool all_valid = fmpz_equal_ui(sum, total_inputs);
        fmpz_clear(sum);
        if (!all_valid) {
            // Batch check fails. Test single.
            // Binary search is more rounds but (possibly) less bytes)
            fmpz_t* sums_other; new_fmpz_array(&sums_other, num_inputs);
            recv_fmpz_batch(serverfd, sums_other, num_inputs);
            for (unsigned int i = 0; i < num_inputs; i++) {
                fmpz_mod_add(sums[i], sums[i], sums_other[i], mod_ctx);
                valid[i] &= fmpz_is_one(sums[i]);
            }
            clear_fmpz_array(sums_other, num_inputs);
        }
        sent_bytes += send_bool_batch(serverfd, valid, num_inputs);

        clear_fmpz_array(sums, num_inputs);
        std::cout << "validate time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        fmpz_t* a; new_fmpz_array(&a, max_inp);
        size_t num_valid = accumulate(num_inputs, max_inp, shares_p, valid, a);
        delete[] valid;
        clear_fmpz_array(shares_p, num_inputs * max_inp);
        std::cout << "accumulate time: " << sec_from(start2) << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;
        start2 = clock_start();

        // receive b
        fmpz_t* b; new_fmpz_array(&b, max_inp);
        recv_fmpz_batch(serverfd, b, max_inp);

        std::cout << "Final valid count: " << num_valid << " / " << total_inputs << std::endl;
        std::cout << "sent server bytes: " << sent_bytes << std::endl;
        if (num_valid < total_inputs * (1 - INVALID_THRESHOLD)) {
            std::cout << "Failing, This is less than the invalid threshold of " << INVALID_THRESHOLD << std::endl;
            clear_fmpz_array(b, max_inp);
            clear_fmpz_array(a, max_inp);
            return RET_INVALID;
        }

        // Sum accumulates
        for (unsigned int j = 0; j < max_inp; j++) {
            fmpz_mod_add(a[j], a[j], b[j], mod_ctx);
        }
        clear_fmpz_array(b, max_inp);
        // output
        for (unsigned int j = 0; j < max_inp && j < 32; j++) {
            std::cout << " Freq(" << j << ") = ";
            fmpz_print(a[j]);
            std::cout << "\n";
        }
        clear_fmpz_array(a, max_inp);
        std::cout << "eval time: " << sec_from(start2) << std::endl;
        return RET_ANS;
    }
}

returnType heavy_op(const initMsg msg, const int clientfd, const int serverfd, const int server_num) {
    auto start = clock_start();

    std::unordered_map<std::string, bool*> share_map;
    int cli_bytes = 0;
    const size_t b = msg.num_bits;
    const unsigned int total_inputs = msg.num_of_inputs;

    // B Not actually necessary for evaluation
    // Could help for "verify" that hashes actually have right sign after finding the output,
    // but it doesn't help evaluation.
    // flint_rand_t hash_seed; flint_randinit(hash_seed);
    // recv_seed(clientfd, hash_seed);
    // HashStore hash_store(b, b, 2, hash_seed);

    fmpz_t* bucket0; new_fmpz_array(&bucket0, b);
    fmpz_t* bucket1; new_fmpz_array(&bucket1, b);

    for (unsigned int i = 0; i < total_inputs; i++) {
        char tag_c[TAG_LENGTH];
        cli_bytes += recv_in(clientfd, &tag_c[0], TAG_LENGTH);
        const std::string tag(tag_c, tag_c+TAG_LENGTH);

        bool* const buff = new bool[2 * b];
        cli_bytes += recv_bool_batch(clientfd, buff, 2 * b);

        if (share_map.find(tag) != share_map.end()) {
            delete[] buff;
            continue;
        }

        share_map[tag] = buff;
    }

    std::cout << "Received " << total_inputs << " total shares" << std::endl;
    std::cout << "bytes from client: " << cli_bytes << std::endl;
    std::cout << "receive time: " << sec_from(start) << std::endl;

    correlated_store->check_DaBits(total_inputs * 2 * b);

    start = clock_start();
    auto start2 = clock_start();
    int sent_bytes = 0;
    size_t num_inputs;

    if (server_num == 1) {
        num_inputs = share_map.size();
        sent_bytes += send_size(serverfd, num_inputs);
    } else {
        recv_size(serverfd, num_inputs);
    }
    const size_t n = num_inputs * b;
    bool* const x = new bool[n];
    bool* const y = new bool[n];
    bool* const valid = new bool[num_inputs];
    memset(valid, true, num_inputs * sizeof(bool));

    /* Share sync */
    if (server_num == 1) {
        size_t idx = 0;
        for (const auto& share : share_map) {
            sent_bytes += send_tag(serverfd, share.first);

            memcpy(&x[idx*b], share.second, b);
            memcpy(&y[idx*b], &(share.second[b]), b);
            idx++;
            delete[] share.second;
        }
        // Sync valid now, since no validation needed
        recv_bool_batch(serverfd, valid, num_inputs);
    } else {
        for (unsigned int i = 0; i < num_inputs; i++) {
            const std::string tag = recv_tag(serverfd);

            valid[i] &= (share_map.find(tag) != share_map.end());
            if (!valid[i]) {
                memset(&x[i * b], 0, b);
                memset(&y[i * b], 0, b);
                continue;
            }
            bool* share = share_map[tag];
            memcpy(&x[i * b], share, b);
            memcpy(&y[i * b], &(share[b]), b);
            delete[] share;
        }
        // Sync valid now, since no validation needed
        sent_bytes += send_bool_batch(serverfd, valid, num_inputs);
    }
    std::cout << "tag time: " << sec_from(start2) << std::endl;
    start2 = clock_start();

    /* Conversion + accumulation */
    sent_bytes += correlated_store->heavy_convert(num_inputs, b, x, y, valid, bucket0, bucket1);

    delete[] x;
    delete[] y;
    std::cout << "convert+accum time: " << sec_from(start2) << std::endl;
    std::cout << "total compute time: " << sec_from(start) << std::endl;
    std::cout << "compute bytes sent: " << sent_bytes << std::endl;
    start2 = clock_start();

    /* Evaluation */
    size_t num_valid = 0;
    for (unsigned int i = 0; i < num_inputs; i++)
        num_valid += valid[i];
    delete[] valid;
    std::cout << "Final valid count: " << num_valid << " / " << total_inputs << std::endl;
    if (num_valid < total_inputs * (1 - INVALID_THRESHOLD)) {
        std::cout << "Failing, This is less than the invalid threshold of " << INVALID_THRESHOLD << std::endl;
        clear_fmpz_array(bucket0, b);
        clear_fmpz_array(bucket1, b);
        return RET_INVALID;
    }
    fmpz_t* larger; new_fmpz_array(&larger, b);

    sent_bytes = correlated_store->abs_cmp(b, bucket0, bucket1, larger);

    clear_fmpz_array(bucket0, b);
    clear_fmpz_array(bucket1, b);

    if (server_num == 1) {
        sent_bytes += send_fmpz_batch(serverfd, larger, b);
        clear_fmpz_array(larger, b);
        std::cout << "evaluate time: " << sec_from(start2) << std::endl;
        std::cout << "evaluate bytes sent: " << sent_bytes << std::endl;

        return RET_NO_ANS;
    } else {
        fmpz_t* larger_other; new_fmpz_array(&larger_other, b);
        recv_fmpz_batch(serverfd, larger_other, b);

        // fmpz_from_bool_array?
        uint64_t ans = 0;
        for (unsigned int j = 0; j < b; j++) {
            fmpz_mod_add(larger[j], larger[j], larger_other[j], mod_ctx);
            const bool bit = fmpz_is_one(larger[j]);
            ans |= ( bit << j );
        }

        std::cout << "### Heavy hitter value is " << ans << std::endl;

        clear_fmpz_array(larger, b);
        clear_fmpz_array(larger_other, b);

        std::cout << "evaluate time: " << sec_from(start2) << std::endl;
        std::cout << "evaluate bytes sent: " << sent_bytes << std::endl;

        return RET_ANS;
    }
}

returnType multi_heavy_op(const initMsg msg, const int clientfd, const int serverfd, const int server_num) {
    auto start = clock_start();

    typedef std::tuple <bool*, bool*, bool*, bool*> sharetype;
    std::unordered_map<std::string, sharetype> share_map;

    int64_t cli_bytes = 0;
    const size_t num_bits = msg.num_bits;
    const unsigned int total_inputs = msg.num_of_inputs;

    // Get other params
    size_t K; recv_size(clientfd, K);
    double delta; recv_double(clientfd, delta);
    double eps; recv_double(clientfd, eps);
    MultiHeavyConfig cfg(K, delta, num_bits, eps, 1);
    cfg.print();

    flint_rand_t hash_seed_classify; flint_randinit(hash_seed_classify);
    flint_rand_t hash_seed_split; flint_randinit(hash_seed_split);
    flint_rand_t hash_seed_value; flint_randinit(hash_seed_value);
    flint_rand_t hash_seed_count; flint_randinit(hash_seed_count);
    recv_seed(clientfd, hash_seed_classify);
    recv_seed(clientfd, hash_seed_split);
    recv_seed(clientfd, hash_seed_value);
    recv_seed(clientfd, hash_seed_count);

    /* Init structures */
    // Q sets of singleHeavy, each with Depth layers. (repeat across B)
    const size_t share_size_sh = cfg.Q * cfg.D;
    // For each Q, identify b (first hash)
    const size_t share_size_mask = cfg.Q * cfg.B;
    // Count-min
    const size_t share_size_count = cfg.countmin_cfg.d * cfg.countmin_cfg.w;
    // SingleHeavy bucket pairs
    const size_t num_sh = cfg.Q * cfg.B * cfg.D;

    if (STORE_TYPE != ot_store) {
        ((CorrelatedStore*) correlated_store)->check_DaBits(
            total_inputs * (share_size_count + share_size_mask + 4 * num_sh));
        ((CorrelatedStore*) correlated_store)->check_BoolTriples(
            3 * total_inputs * num_sh);
    }

    // Get client
    auto start2 = clock_start();
    for (unsigned int i = 0; i < total_inputs; i++) {
        char tag_c[TAG_LENGTH];
        cli_bytes += recv_in(clientfd, &tag_c[0], TAG_LENGTH);
        const std::string tag(tag_c, tag_c+TAG_LENGTH);

        bool* const shares_sh_x = new bool[share_size_sh];
        bool* const shares_sh_y = new bool[share_size_sh];
        bool* const shares_mask = new bool[share_size_mask];
        bool* const shares_count = new bool[share_size_count];
        cli_bytes += recv_bool_batch(clientfd, shares_sh_x, share_size_sh);
        cli_bytes += recv_bool_batch(clientfd, shares_sh_y, share_size_sh);
        cli_bytes += recv_bool_batch(clientfd, shares_mask, share_size_mask);
        cli_bytes += recv_bool_batch(clientfd, shares_count, share_size_count);

        if (share_map.find(tag) != share_map.end()) {
            delete[] shares_sh_x;
            delete[] shares_sh_y;
            delete[] shares_mask;
            delete[] shares_count;
            continue;
        }

        share_map[tag] = {shares_sh_x, shares_sh_y, shares_mask, shares_count};
    }

    std::cout << "Received " << total_inputs << " total shares" << std::endl;
    std::cout << "bytes from client: " << cli_bytes << std::endl;
    std::cout << "receive time: " << sec_from(start2) << std::endl;

    start = clock_start();
    start2 = clock_start();
    int64_t sent_bytes = 0;

    /* Sync shares */
    size_t num_inputs;
    if (server_num == 1) {
        num_inputs = share_map.size();
        sent_bytes += send_size(serverfd, num_inputs);
    } else {
        recv_size(serverfd, num_inputs);
    }

    std::cout << "num_inputs: " << num_inputs << std::endl;
    bool* const shares_sh_x = new bool[num_inputs * share_size_sh];
    bool* const shares_sh_y = new bool[num_inputs * share_size_sh];
    bool* const shares_mask = new bool[num_inputs * share_size_mask];
    bool* const shares_count = new bool[num_inputs * share_size_count];
    bool* const valid = new bool[num_inputs];
    memset(valid, true, num_inputs * sizeof(bool));

    if (server_num == 1) {
        size_t i = 0;
        for (const auto& share : share_map) {
            sent_bytes += send_tag(serverfd, share.first);

            bool* a; bool* b; bool* c; bool* d;
            std::tie(a, b, c, d) = share.second;
            memcpy(&shares_sh_x[i * share_size_sh], a, share_size_sh);
            memcpy(&shares_sh_y[i * share_size_sh], b, share_size_sh);
            memcpy(&shares_mask[i * share_size_mask], c, share_size_mask);
            memcpy(&shares_count[i * share_size_count], d, share_size_count);
            delete a;
            delete b;
            delete c;
            delete d;

            i++;
        }
    } else {
        for (unsigned int i = 0; i < num_inputs; i++) {
            const std::string tag = recv_tag(serverfd);

            valid[i] &= (share_map.find(tag) != share_map.end());
            if (!valid[i]) continue;

            bool* a; bool* b; bool* c; bool* d;
            std::tie(a, b, c, d) = share_map[tag];
            memcpy(&shares_sh_x[i * share_size_sh], a, share_size_sh);
            memcpy(&shares_sh_y[i * share_size_sh], b, share_size_sh);
            memcpy(&shares_mask[i * share_size_mask], c, share_size_mask);
            memcpy(&shares_count[i * share_size_count], d, share_size_count);
            delete a;
            delete b;
            delete c;
            delete d;
        }
    }
    std::cout << "tag time: " << sec_from(start2) << std::endl;
    start2 = clock_start();
    auto start3 = clock_start();

    /* Conversion, accumulation */
    /* Stages:
    Validate each b (shares_mask) is frequency vector. With B2A to check sum to 1
    "Validate H_q vector correct". shares_sh. Not needed for R=1, is just entry
    (joint share validity)
    Update Countmin with shares_count (and freq check in parallel)
    For each (q, b):
      final mask = shares_sh & R_mask.
      SH update, except OT multiply [z] by mask.

    Round 1: B2A count (use+valid), mask (valid), first convert mult
    Round 2: freq valid values, second convert mult
    Round 3: freq valid results, convert B2A

    mask: B2A just for validation. Kept as bool for Ot multiply
    sh: half B2A for use (z), half kept selection
    count: B2A for validation and accumulation
    */

    // 2 * num_inputs, 4*len_all + num*(count + mask)
    const size_t len_all = num_inputs * cfg.Q * cfg.B * cfg.D;
    const size_t len_p = num_inputs * (share_size_count + share_size_mask);
    bool* const send_buff = new bool[4 * len_all + len_p];
    bool* const recv_buff = new bool[4 * len_all + len_p];
    memset(recv_buff, 0, 4 * len_all + len_p);

    // Setup convert stage 1
    bool* const m_ext = new bool[2*len_all];
    bool* const xy_ext = new bool[2*len_all];
    cross_fill_bool(num_inputs * cfg.Q, cfg.B, cfg.D,
            shares_mask, shares_sh_x, m_ext, xy_ext);
    delete[] shares_sh_x;
    cross_fill_bool(num_inputs * cfg.Q, cfg.B, cfg.D,
            shares_mask, shares_sh_y, &m_ext[len_all], &xy_ext[len_all]);
    delete[] shares_sh_y;
    bool* const z = new bool[4*len_all];
    memcpy(z, m_ext, len_all);
    // 4*len_all into send buff
    correlated_store->multiply_BoolShares_setup(2*len_all, xy_ext, m_ext,
            &z[len_all], send_buff);
    // B2A setup
    fmpz_t* shares_p; new_fmpz_array(&shares_p, len_p);
    correlated_store->b2a_single_setup(num_inputs * share_size_count,
            shares_count, shares_p, &send_buff[4*len_all]);
    delete[] shares_count;
    correlated_store->b2a_single_setup(num_inputs * share_size_mask,
            shares_mask, &shares_p[num_inputs * share_size_count],
            &send_buff[4*len_all + num_inputs * share_size_count]);
    delete[] shares_mask;

    /* Round 1: convert mult 1, B2A count, B2A mask */
    sent_bytes += swap_bool_batch(serverfd, send_buff, recv_buff,
            4*len_all + len_p);
    std::cout << "Round 1 time: " << sec_from(start3) << std::endl;
    start3 = clock_start();

    // Finish B2A
    correlated_store->b2a_single_finish(num_inputs * share_size_count,
            shares_p, &send_buff[4*len_all], &recv_buff[4*len_all]);
    correlated_store->b2a_single_finish(num_inputs * share_size_mask,
            &shares_p[num_inputs * share_size_count],
            &send_buff[4*len_all + num_inputs * share_size_count],
            &recv_buff[4*len_all + num_inputs * share_size_count]);
    // Finish convert stage 1, setup convert stage 2
    correlated_store->multiply_BoolShares_finish(2*len_all, xy_ext, m_ext,
            &z[len_all], send_buff, recv_buff);
    delete[] m_ext;
    correlated_store->multiply_BoolShares_setup(len_all, xy_ext,
            &z[2*len_all], &z[3*len_all], send_buff);
    fmpz_t* sums; new_fmpz_array(&sums, cfg.Q + cfg.countmin_cfg.d);
    sum_accum(num_inputs, cfg.countmin_cfg.d, cfg.countmin_cfg.w, shares_p, sums);
    sum_accum(num_inputs, cfg.Q, cfg.B,
            &shares_p[num_inputs * share_size_count], &sums[cfg.countmin_cfg.d]);
    fmpz_t* sums_other; new_fmpz_array(&sums_other, cfg.Q + cfg.countmin_cfg.d);

    /* Round 2: freq valid values, second convert mult */
    sent_bytes += swap_bool_fmpz_batch(serverfd,
            send_buff, recv_buff, 2 * len_all,
            sums, sums_other, cfg.Q + cfg.countmin_cfg.d);
    std::cout << "Round 2 time: " << sec_from(start3) << std::endl;
    start3 = clock_start();

    // Finish convert stage 2
    correlated_store->multiply_BoolShares_finish(len_all, xy_ext,
            &z[2*len_all], &z[3*len_all], send_buff, recv_buff);
    delete[] xy_ext;
    // Finish freq, update valid
    bool all_valid = true;
    for (unsigned int j = 0; j < cfg.Q + cfg.countmin_cfg.d; j++) {
        fmpz_mod_add(sums[j], sums[j], sums_other[j], mod_ctx);
        all_valid &= fmpz_equal_ui(sums[j], num_inputs);
    }
    clear_fmpz_array(sums, cfg.Q + cfg.countmin_cfg.d);
    clear_fmpz_array(sums_other, cfg.Q + cfg.countmin_cfg.d);
    if (!all_valid) {
        memset(valid, false, num_inputs);
        std::cout << "Batch not valid. Individual check currently not implemented\n";
    }
    memcpy(send_buff, valid, num_inputs);
    // setup convert stage 3
    fmpz_t* zp; new_fmpz_array(&zp, 4 * len_all);
    std::cout << "b2a setup 2: " << 4 * len_all << std::endl;
    correlated_store->b2a_single_setup(4 * len_all, z, zp, &send_buff[num_inputs]);
    delete[] z;

    /* Round 3: convert stage 3 (b2a), valid swap */
    sent_bytes += swap_bool_batch(serverfd, send_buff, recv_buff,
            num_inputs + 4 * len_all);
    std::cout << "Round 3 time: " << sec_from(start3) << std::endl;
    start3 = clock_start();

    // Finish heavy convert
    correlated_store->b2a_single_finish(4 * len_all, zp,
            &send_buff[num_inputs], &recv_buff[num_inputs]);
    delete[] send_buff;
    // Merge valid
    for (unsigned int i = 0; i < num_inputs; i++) {
        valid[i] &= recv_buff[i];
    }
    delete[] recv_buff;
    // Heavy accum
    fmpz_t* bucket0; new_fmpz_array(&bucket0, num_sh);
    fmpz_t* bucket1; new_fmpz_array(&bucket1, num_sh);
    correlated_store->heavy_accumulate(num_inputs, cfg.Q * cfg.B * cfg.D, zp,
            valid, bucket0, bucket1);
    clear_fmpz_array(zp, 4 * len_all);
    // Count-min accum
    fmpz_t* countmin_accum; new_fmpz_array(&countmin_accum, share_size_count);
    accumulate(num_inputs, share_size_count, shares_p, valid, countmin_accum);
    clear_fmpz_array(shares_p, len_p);

    std::cout << "accum time: " << sec_from(start3) << std::endl;
    std::cout << "total accum time: " << sec_from(start2) << std::endl;
    std::cout << "total compute time: " << sec_from(start) << std::endl;
    std::cout << "compute bytes sent: " << sent_bytes << std::endl;
    start2 = clock_start();

    /* Evaluation */
    size_t num_valid = 0;
    for (unsigned int i = 0; i < num_inputs; i++)
        num_valid += valid[i];
    std::cout << "Final valid count: " << num_valid << " / " << total_inputs << std::endl;
    delete[] valid;
    if (num_valid < total_inputs * (1 - INVALID_THRESHOLD)) {
        std::cout << "Failing, This is less than the invalid threshold of " << INVALID_THRESHOLD << std::endl;
        clear_fmpz_array(bucket0, num_sh);
        clear_fmpz_array(bucket1, num_sh);
        clear_fmpz_array(countmin_accum, share_size_count);
        return RET_INVALID;
    }

    uint64_t* const top_values = new uint64_t[K];
    uint64_t* const top_freqs = new uint64_t[K];

    full_heavy_extract(server_num, cfg, bucket0, bucket1, hash_seed_split, hash_seed_count,
            countmin_accum, num_inputs, top_values, top_freqs);
    garbleIO->flush();

    clear_fmpz_array(bucket0, num_sh);
    clear_fmpz_array(bucket1, num_sh);
    // Not needed, since freed within extract
    // clear_fmpz_array(countmin_accum, share_size_count);

    std::cout << "Top K = " << K << " values and freqs, decreasing\n";
    for (unsigned int i = 0; i < K; i++) {
      std::cout << "Value: " << top_values[i] << ", freq: " << top_freqs[i] << std::endl;
    }
    std::cout << "evaluate time: " << sec_from(start2) << std::endl;

    delete[] top_values;
    delete[] top_freqs;

    return RET_ANS;
}

returnType top_k_op(const initMsg msg, const int clientfd, const int serverfd, const int server_num) {
    auto start = clock_start();

    typedef std::tuple <bool*, bool*, bool*, bool*, bool*> sharetype;
    std::unordered_map<std::string, sharetype> share_map;

    int64_t cli_bytes = 0;
    const size_t num_bits = msg.num_bits;
    const unsigned int total_inputs = msg.num_of_inputs;

    // Get other params
    size_t K; recv_size(clientfd, K);
    double delta; recv_double(clientfd, delta);
    double eps; recv_double(clientfd, eps);
    size_t R; recv_size(clientfd, R);
    MultiHeavyConfig cfg(K, delta, num_bits, eps, R);
    cfg.print();

    flint_rand_t hash_seed_classify; flint_randinit(hash_seed_classify);
    flint_rand_t hash_seed_split; flint_randinit(hash_seed_split);
    flint_rand_t hash_seed_value; flint_randinit(hash_seed_value);
    flint_rand_t hash_seed_count; flint_randinit(hash_seed_count);
    recv_seed(clientfd, hash_seed_classify);
    recv_seed(clientfd, hash_seed_split);
    recv_seed(clientfd, hash_seed_value);
    recv_seed(clientfd, hash_seed_count);
    // R seed?

    /* Init structures */
    // Q sets of singleHeavy, each with Depth layers. (repeat across B)
    const size_t share_size_sh = cfg.Q * cfg.D;
    // For each Q, identify b (first hash)
    const size_t share_size_bucket = cfg.Q * cfg.B;
    // Layer selector.
    const size_t share_size_layer = cfg.R;
    // Count-min
    const size_t share_size_count = cfg.countmin_cfg.d * cfg.countmin_cfg.w;
    // SingleHeavy bucket pairs
    const size_t num_sh = cfg.R * cfg.Q * cfg.B * cfg.D;

    if (STORE_TYPE != ot_store) {
        ((CorrelatedStore*) correlated_store)->check_DaBits(
                total_inputs * (share_size_count + share_size_bucket + 4 * num_sh));
        ((CorrelatedStore*) correlated_store)->check_BoolTriples(
                total_inputs * (3 * num_sh + share_size_sh
                                + share_size_bucket * share_size_layer));
    }

    fmpz_t* bucket0; new_fmpz_array(&bucket0, num_sh);
    fmpz_t* bucket1; new_fmpz_array(&bucket1, num_sh);

    // Get client
    auto start2 = clock_start();
    for (unsigned int i = 0; i < total_inputs; i++) {
        char tag_c[TAG_LENGTH];
        cli_bytes += recv_in(clientfd, &tag_c[0], TAG_LENGTH);
        const std::string tag(tag_c, tag_c+TAG_LENGTH);

        bool* const shares_sh_x = new bool[share_size_sh];
        bool* const shares_sh_y = new bool[share_size_sh];
        bool* const shares_bucket = new bool[share_size_bucket];
        bool* const shares_layer = new bool[share_size_layer];
        bool* const shares_count = new bool[share_size_count];
        cli_bytes += recv_bool_batch(clientfd, shares_sh_x, share_size_sh);
        cli_bytes += recv_bool_batch(clientfd, shares_sh_y, share_size_sh);
        cli_bytes += recv_bool_batch(clientfd, shares_bucket, share_size_bucket);
        cli_bytes += recv_bool_batch(clientfd, shares_layer, share_size_layer);
        cli_bytes += recv_bool_batch(clientfd, shares_count, share_size_count);

        if (share_map.find(tag) != share_map.end()) {
            delete[] shares_sh_x;
            delete[] shares_sh_y;
            delete[] shares_bucket;
            delete[] shares_layer;
            delete[] shares_count;
            continue;
        }

        share_map[tag] = {shares_sh_x, shares_sh_y, shares_bucket, shares_layer, shares_count};
    }

    std::cout << "Received " << total_inputs << " total shares" << std::endl;
    std::cout << "bytes from client: " << cli_bytes << std::endl;
    std::cout << "receive time: " << sec_from(start2) << std::endl;

    start = clock_start();
    start2 = clock_start();
    int64_t sent_bytes = 0;

    /* Sync shares */
    size_t num_inputs;
    if (server_num == 1) {
        num_inputs = share_map.size();
        sent_bytes += send_size(serverfd, num_inputs);
    } else {
        recv_size(serverfd, num_inputs);
    }

    std::cout << "num_inputs: " << num_inputs << std::endl;
    bool* const shares_sh_x = new bool[num_inputs * share_size_sh];
    bool* const shares_sh_y = new bool[num_inputs * share_size_sh];
    bool* const shares_bucket = new bool[num_inputs * share_size_bucket];
    bool* const shares_layer = new bool[num_inputs * share_size_layer];
    bool* const shares_count = new bool[num_inputs * share_size_count];
    bool* const valid = new bool[num_inputs];
    memset(valid, true, num_inputs * sizeof(bool));

    if (server_num == 1) {
        size_t i = 0;
        for (const auto& share : share_map) {
            sent_bytes += send_tag(serverfd, share.first);

            bool* a; bool* b; bool* c; bool* d; bool* e;
            std::tie(a, b, c, d, e) = share.second;
            memcpy(&shares_sh_x[i * share_size_sh], a, share_size_sh);
            memcpy(&shares_sh_y[i * share_size_sh], b, share_size_sh);
            memcpy(&shares_bucket[i * share_size_bucket], c, share_size_bucket);
            memcpy(&shares_layer[i * share_size_layer], d, share_size_layer);
            memcpy(&shares_count[i * share_size_count], e, share_size_count);
            delete a;
            delete b;
            delete c;
            delete d;
            delete e;

            i++;
        }
    } else {
        for (unsigned int i = 0; i < num_inputs; i++) {
            const std::string tag = recv_tag(serverfd);

            valid[i] &= (share_map.find(tag) != share_map.end());
            if (!valid[i]) continue;

            bool* a; bool* b; bool* c; bool* d; bool* e;
            std::tie(a, b, c, d, e) = share_map[tag];
            memcpy(&shares_sh_x[i * share_size_sh], a, share_size_sh);
            memcpy(&shares_sh_y[i * share_size_sh], b, share_size_sh);
            memcpy(&shares_bucket[i * share_size_bucket], c, share_size_bucket);
            memcpy(&shares_layer[i * share_size_layer], d, share_size_layer);
            memcpy(&shares_count[i * share_size_count], e, share_size_count);
            delete a;
            delete b;
            delete c;
            delete d;
            delete e;
        }
    }
    std::cout << "tag time: " << sec_from(start2) << std::endl;
    start2 = clock_start();
    auto start3 = clock_start();

    /*
    Round 1: B2A count (use+valid), mask (valid), first convert mult
    Round 2: freq valid values, second convert mult
    Round 3: freq valid results, convert B2A

    mask: B2A just for validation. Kept as bool for Ot multiply
    sh: half B2A for use (z), half kept selection
    count: B2A for validation and accumulation
    */

    // Length declares
    const size_t len_all = num_inputs * cfg.Q * cfg.B * cfg.R * cfg.D;
    const size_t len_part = num_inputs * cfg.Q * (cfg.D + cfg.B * cfg.R);
    const size_t len_p = num_inputs * (share_size_count + share_size_bucket);
    bool* const send_buff = new bool[6 * len_all];
    bool* const recv_buff = new bool[6 * len_all];
    memset(recv_buff, 0, 6 * len_all);

    // Setup convert stage 1: xy normal, cross B, R
    bool* const a = new bool[3 * len_all];
    bool* const b = new bool[3 * len_all];
    bool* const c = new bool[len_part];
    memcpy(a, shares_sh_x, num_inputs * cfg.Q * cfg.D);
    memcpy(b, shares_sh_y, num_inputs * cfg.Q * cfg.D);
    cross_fill_bool(num_inputs * cfg.Q, cfg.B, cfg.R, shares_bucket, shares_layer,
            &a[num_inputs * cfg.Q * cfg.D], &b[num_inputs * cfg.Q * cfg.D]);
    delete[] shares_layer;
    correlated_store->multiply_BoolShares_setup(len_part, a, b, c, send_buff);
    // 2 * len_part into send buff
    fmpz_t* shares_p; new_fmpz_array(&shares_p, len_p);
    correlated_store->b2a_single_setup(num_inputs * share_size_count,
            shares_count, shares_p, &send_buff[2*len_part]);
    delete[] shares_count;
    correlated_store->b2a_single_setup(num_inputs * share_size_bucket,
            shares_bucket, &shares_p[num_inputs * share_size_count],
            &send_buff[2*len_part + num_inputs * share_size_count]);
    delete[] shares_bucket;

    /* Round 1: convert mult 1, B2A count, B2A mask */
    sent_bytes += swap_bool_batch(serverfd, send_buff, recv_buff,
            2*len_part + len_p);
    std::cout << "Round 1 time: " << sec_from(start3) << std::endl;
    start3 = clock_start();

    // Finish B2A
    correlated_store->b2a_single_finish(num_inputs * share_size_count,
            shares_p, &send_buff[2*len_part], &recv_buff[2*len_part]);
    correlated_store->b2a_single_finish(num_inputs * share_size_bucket,
            &shares_p[num_inputs * share_size_count],
            &send_buff[2*len_part + num_inputs * share_size_count],
            &recv_buff[2*len_part + num_inputs * share_size_count]);
    // Finish convert stage 1, setup convert stage 2
    correlated_store->multiply_BoolShares_finish(len_part, a, b, c, send_buff, recv_buff);
    // (x, y, xy) * m1m2
    cross_fill_bool(num_inputs * cfg.Q, cfg.B * cfg.R, cfg.D,
            &c[num_inputs * cfg.Q * cfg.D], shares_sh_x, b, a);
    delete[] shares_sh_x;
    cross_fill_bool(num_inputs * cfg.Q, cfg.B * cfg.R, cfg.D,
            &c[num_inputs * cfg.Q * cfg.D], shares_sh_y, &b[len_all], &a[len_all]);
    delete[] shares_sh_y;
    cross_fill_bool(num_inputs * cfg.Q, cfg.B * cfg.R, cfg.D,
            &c[num_inputs * cfg.Q * cfg.D], c, &b[2*len_all], &a[2*len_all]);
    delete[] c;
    bool* const z = new bool[4 * len_all];
    memcpy(z, b, len_all);
    correlated_store->multiply_BoolShares_setup(3 * len_all, a, b, &z[len_all], send_buff);
    fmpz_t* sums; new_fmpz_array(&sums, cfg.Q + cfg.countmin_cfg.d);
    sum_accum(num_inputs, cfg.countmin_cfg.d, cfg.countmin_cfg.w, shares_p, sums);
    sum_accum(num_inputs, cfg.Q, cfg.B,
            &shares_p[num_inputs * share_size_count], &sums[cfg.countmin_cfg.d]);
    fmpz_t* sums_other; new_fmpz_array(&sums_other, cfg.Q + cfg.countmin_cfg.d);

    /* Round 2: Freq valid values, second convert mult */
    sent_bytes += swap_bool_fmpz_batch(serverfd,
            send_buff, recv_buff, 2 * 3 * len_all,
            sums, sums_other, cfg.Q + cfg.countmin_cfg.d);
    std::cout << "Round 2 time: " << sec_from(start3) << std::endl;
    start3 = clock_start();

    // Finish convert stage 2
    correlated_store->multiply_BoolShares_finish(3 * len_all, a, b, &z[len_all], send_buff, recv_buff);
    delete[] a;
    delete[] b;
    // Finish freq, update valid
    bool all_valid = true;
    for (unsigned int j = 0; j < cfg.Q + cfg.countmin_cfg.d; j++) {
        fmpz_mod_add(sums[j], sums[j], sums_other[j], mod_ctx);
        all_valid &= fmpz_equal_ui(sums[j], num_inputs);
    }
    clear_fmpz_array(sums, cfg.Q + cfg.countmin_cfg.d);
    clear_fmpz_array(sums_other, cfg.Q + cfg.countmin_cfg.d);
    if (!all_valid) {
        memset(valid, false, num_inputs);
        std::cout << "Batch not valid. Individual check currently not implemented\n";
    }
    memcpy(send_buff, valid, num_inputs);
    // setup convert stage 3
    fmpz_t* zp; new_fmpz_array(&zp, 4 * len_all);
    correlated_store->b2a_single_setup(4 * len_all, z, zp, &send_buff[num_inputs]);
    delete[] z;

    /* Round 3: convert stage 3, valid swap */
    sent_bytes += swap_bool_batch(serverfd, send_buff, recv_buff,
            num_inputs + 4 * len_all);
    std::cout << "Round 3 time: " << sec_from(start3) << std::endl;
    start3 = clock_start();

    // Finish heavy convert
    correlated_store->b2a_single_finish(4 * len_all, zp, &send_buff[num_inputs], &recv_buff[num_inputs]);
    delete[] send_buff;
    // Merge valid
    for (unsigned int i = 0; i < num_inputs; i++) {
        valid[i] &= recv_buff[i];
    }
    delete[] recv_buff;
    // Heavy accum
    correlated_store->heavy_accumulate(num_inputs,
            cfg.Q * cfg.B * cfg.R * cfg.D, zp, valid, bucket0, bucket1);
    clear_fmpz_array(zp, 4 * len_all);
    // Count-min accum
    fmpz_t* countmin_accum; new_fmpz_array(&countmin_accum, share_size_count);
    accumulate(num_inputs, share_size_count, shares_p, valid, countmin_accum);
    clear_fmpz_array(shares_p, len_p);

    std::cout << "accum time: " << sec_from(start3) << std::endl;
    std::cout << "total accum time: " << sec_from(start2) << std::endl;
    std::cout << "total compute time: " << sec_from(start) << std::endl;
    std::cout << "compute bytes sent: " << sent_bytes << std::endl;
    start2 = clock_start();

    /* Evaluation */
    size_t num_valid = 0;
    for (unsigned int i = 0; i < num_inputs; i++)
        num_valid += valid[i];
    std::cout << "Final valid count: " << num_valid << " / " << total_inputs << std::endl;
    delete[] valid;
    if (num_valid < total_inputs * (1 - INVALID_THRESHOLD)) {
        std::cout << "Failing, This is less than the invalid threshold of " << INVALID_THRESHOLD << std::endl;
        clear_fmpz_array(bucket0, num_sh);
        clear_fmpz_array(bucket1, num_sh);
        clear_fmpz_array(countmin_accum, share_size_count);
        return RET_INVALID;
    }

    uint64_t* const top_values = new uint64_t[K];
    uint64_t* const top_freqs = new uint64_t[K];

    full_heavy_extract(server_num, cfg, bucket0, bucket1, hash_seed_split, hash_seed_count,
            countmin_accum, num_inputs, top_values, top_freqs);
    garbleIO->flush();

    clear_fmpz_array(bucket0, num_sh);
    clear_fmpz_array(bucket1, num_sh);
    // Not needed, since freed within extract
    // clear_fmpz_array(countmin_accum, share_size_count);

    std::cout << "Top K = " << K << " values and freqs, decreasing\n";
    for (unsigned int i = 0; i < K; i++) {
      std::cout << "Value: " << top_values[i] << ", freq: " << top_freqs[i] << std::endl;
    }
    std::cout << "evaluate time: " << sec_from(start2) << std::endl;

    delete[] top_values;
    delete[] top_freqs;

    return RET_ANS;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cout << "Usage: ./bin/server server_num(0/1) this_client_port server0_port" << endl;
        return 1;
    }

    const int server_num = atoi(argv[1]);  // Server # 1 or # 2
    const int client_port = atoi(argv[2]); // port of this server, for the client
    const int server_port = atoi(argv[3]); // port of this server, for the other server

    std::cout << "This server is server # " << server_num << std::endl;
    std::cout << "  Listening for client on " << client_port << std::endl;
    std::cout << "  Listening for server on " << server_port << std::endl;

    init_constants();

    // Set serverfd
    // Server 0 listens, via newsockfd_server
    // Server 1 connects, via sockfd_server
    int sockfd_server, newsockfd_server, serverfd = 0;
    if (server_num == 0) {
        server0_listen(sockfd_server, newsockfd_server, server_port, 1);
        serverfd = newsockfd_server;
    } else if (server_num == 1) {
        server1_connect(sockfd_server, server_port, 1);
        serverfd = sockfd_server;
    } else {
        error_exit("Can only handle servers #0 and #1");
    }

    mult_eval_manager = new MultEvalManager(server_num, serverfd);

    syncSnipSeeds(serverfd, server_num);

    ot0 = new OT_Wrapper(server_num == 0 ? nullptr : SERVER0_IP, 60051);
    ot1 = new OT_Wrapper(server_num == 1 ? nullptr : SERVER1_IP, 60052);

    // TODO: OT disabled/not supported for now.
    if (STORE_TYPE == precompute_store) {
        correlated_store = new PrecomputeStore(
                serverfd, server_num, ot0, ot1, CACHE_SIZE, lazy_precompute);
    // } else if (STORE_TYPE == ot_store) {
    //     correlated_store = new OTCorrelatedStore(serverfd, server_num, ot0, ot1);
    } else if (STORE_TYPE == validate_store) {
        correlated_store = new ValidateCorrelatedStore(
                serverfd, server_num, ot0, ot1, CACHE_SIZE, lazy_precompute);
    } else {
        error_exit("Unknown/unsupported store type");
    }

    // Reuse IO from OT
    garbleIO = ot0->io;
    setup_semi_honest(garbleIO, server_num + 1);
    garbleIO->flush();

    int sockfd, newsockfd;
    sockaddr_in addr;


    bind_and_listen(addr, sockfd, client_port, 1);

    while(1) {
        if (STORE_TYPE != ot_store) {
            ((PrecomputeStore*) correlated_store)->maybe_Update();
        }
        correlated_store->print_Sizes();

        socklen_t addrlen = sizeof(addr);

        std::cout << "waiting for connection..." << std::endl;

        newsockfd = accept(sockfd, (struct sockaddr*)&addr, &addrlen);
        if (newsockfd < 0) error_exit("Connection creation failure");

        // Get an initMsg
        initMsg msg;
        recv_in(newsockfd, &msg, sizeof(initMsg));

        if (msg.type == BIT_SUM_OP) {
            std::cout << "BIT_SUM" << std::endl;
            auto start = clock_start();

            uint64_t ans;
            returnType ret = bit_sum(msg, newsockfd, serverfd, server_num, ans);
            if (ret == RET_ANS)
                std::cout << "Ans: " << ans << std::endl;

            std::cout << "Total time  : " << sec_from(start) << std::endl;
        } else if (msg.type == INT_SUM_OP) {
            std::cout << "INT_SUM" << std::endl;
            auto start = clock_start();

            uint64_t ans;
            returnType ret = int_sum(msg, newsockfd, serverfd, server_num, ans);
            if (ret == RET_ANS)
                std::cout << "Ans: " << ans << std::endl;

            std::cout << "Total time  : " << sec_from(start) << std::endl;
        } else if (msg.type == AND_OP) {
            std::cout << "AND_OP" << std::endl;
            auto start = clock_start();

            bool ans;
            returnType ret = xor_op(msg, newsockfd, serverfd, server_num, ans);
            if (ret == RET_ANS)
                std::cout << "Ans: " << std::boolalpha << ans << std::endl;

            std::cout << "Total time  : " << sec_from(start) << std::endl;
        } else if (msg.type == OR_OP) {
            std::cout << "OR_OP" << std::endl;
            auto start = clock_start();

            bool ans;
            returnType ret = xor_op(msg, newsockfd, serverfd, server_num, ans);
            if (ret == RET_ANS)
                std::cout << "Ans: " << std::boolalpha << ans << std::endl;

            std::cout << "Total time  : " << sec_from(start) << std::endl;
        } else if (msg.type == MAX_OP) {
            std::cout << "MAX_OP" << std::endl;
            auto start = clock_start();

            uint64_t ans;
            returnType ret = max_op(msg, newsockfd, serverfd, server_num, ans);
            if (ret == RET_ANS)
                std::cout << "Ans: " << ans << std::endl;

            std::cout << "Total time  : " << sec_from(start) << std::endl;
        } else if (msg.type == MIN_OP) {
            std::cout << "MIN_OP" << std::endl;
            auto start = clock_start();

            uint64_t ans;
            returnType ret = max_op(msg, newsockfd, serverfd, server_num, ans);
            if (ret == RET_ANS)
                std::cout << "Ans: " << ans << std::endl;

            std::cout << "Total time  : " << sec_from(start) << std::endl;
        } else if (msg.type == VAR_OP) {
            std::cout << "VAR_OP" << std::endl;
            auto start = clock_start();

            double ans;
            returnType ret = var_op(msg, newsockfd, serverfd, server_num, ans);
            if (ret == RET_ANS)
                std::cout << "Ans: " << ans << std::endl;

            std::cout << "Total time  : " << sec_from(start) << std::endl;
        } else if (msg.type == STDDEV_OP) {
            std::cout << "STDDEV_OP" << std::endl;
            auto start = clock_start();

            double ans;
            returnType ret = var_op(msg, newsockfd, serverfd, server_num, ans);
            if (ret == RET_ANS)
                std::cout << "Ans: " << ans << std::endl;

            std::cout << "Total time  : " << sec_from(start) << std::endl;
        } else if (msg.type == LINREG_OP) {
            std::cout << "LINREG_OP" << std::endl;
            auto start = clock_start();

            returnType ret = linreg_op(msg, newsockfd, serverfd, server_num);
            if (ret == RET_ANS)
                ;  // Answer output by linreg_op

            std::cout << "Total time  : " << sec_from(start) << std::endl;
        } else if (msg.type == FREQ_OP) {
            std::cout << "FREQ_OP" << std::endl;
            auto start = clock_start();

            returnType ret = freq_op(msg, newsockfd, serverfd, server_num);
            if (ret == RET_ANS)
                ; // Answer output by freq_op

            std::cout << "Total time  : " << sec_from(start) << std::endl;
        } else if (msg.type == HEAVY_OP) {
            std::cout << "HEAVY_OP" << std::endl;
            auto start = clock_start();

            returnType ret = heavy_op(msg, newsockfd, serverfd, server_num);
            if (ret == RET_ANS)
                ; // Answer output by heavy_op

            std::cout << "Total time  : " << sec_from(start) << std::endl;
        } else if (msg.type == MULTI_HEAVY_OP) {
            std::cout << "MULTI_HEAVY_OP" << std::endl;
            auto start = clock_start();

            returnType ret = multi_heavy_op(msg, newsockfd, serverfd, server_num);
            if (ret == RET_ANS)
                ; // Answer output by multi_heavy_op

            std::cout << "Total time  : " << sec_from(start) << std::endl;
        } else if (msg.type == TOP_K_OP) {
            std::cout << "TOPK_OP" << std::endl;
            auto start = clock_start();

            returnType ret = top_k_op(msg, newsockfd, serverfd, server_num);
            if (ret == RET_ANS)
                ; // Answer output by top_k_op

            std::cout << "Total time  : " << sec_from(start) << std::endl;
        } else if (msg.type == NONE_OP) {
            std::cout << "Empty client message" << std::endl;
        } else {
            std::cout << "Unrecognized message type: " << msg.type << std::endl;
        }
        close(newsockfd);
    }

    delete correlated_store;

    delete ot0;
    delete ot1;

    delete mult_eval_manager;
    RootManager(1).clearCache();
    finalize_semi_honest();
    clear_constants();

    return 0;
}
