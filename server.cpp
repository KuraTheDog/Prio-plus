#include "server.h"

#include <netinet/in.h>
#include <sys/socket.h>
#include <sys/wait.h>
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
#define LAZY_PRECOMPUTE false

// Note: Currently does it in a single batch.
// I.e. recieve and store all, then process all.
// TODO: Set up batching. Either every N inputs (based on space consumed), or every S seconds (figure out timer)

size_t send_out(const int sockfd, const void* const buf, const size_t len) {
    size_t ret = send(sockfd, buf, len, 0);
    if (ret <= 0) error_exit("Failed to send");
    return ret;
}

size_t send_tag(const int sockfd, const std::string x) {
    return send_out(sockfd, &x[0], TAG_LENGTH);
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

std::string get_tag(const int serverfd) {
    char tag_buf[TAG_LENGTH];
    recv_in(serverfd, &tag_buf[0], TAG_LENGTH);
    std::string tag(tag_buf, tag_buf + TAG_LENGTH);
    return tag;
}

int recv_unvalidated(const int clientfd, const std::string tag, const size_t n) {
    if (STORE_TYPE != validate_store) {
        return 0;
    }

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
    if (STORE_TYPE != validate_store) {
        return;
    }
    ((ValidateCorrelatedStore*) correlated_store)->process_Unvalidated(tag);
}

// B2A Multi with dimensions, and int share_2s
// Currently shares_2 and shares_p are flat num_shares*num_values array.
// TODO: split, for the sake of communication
// TODO: move to store?
int const share_convert(const size_t num_shares,  // # inputs
                        const size_t num_values,  // # values per input
                        const size_t* const num_bits,  // # bits per value
                        const uint64_t* const shares_2,
                        fmpz_t* const shares_p
                        ) {
    auto start = clock_start();
    int sent_bytes = 0;

    // convert
    fmpz_t* f_shares2; new_fmpz_array(&f_shares2, num_shares * num_values);
    for (unsigned int i = 0; i < num_shares; i++) {
        for (unsigned int j = 0; j < num_values; j++)
            fmpz_set_ui(f_shares2[i * num_values + j],
                        shares_2[i * num_values + j]);
    }

    size_t* const bits_arr = new size_t[num_shares * num_values];
    // num_bits 1 char, so this is fine
    for (unsigned int i = 0; i < num_shares; i++)
        memcpy(&bits_arr[i * num_values], num_bits, num_values * sizeof(size_t));

    sent_bytes += correlated_store->b2a_multi(
        num_shares * num_values, bits_arr, f_shares2, shares_p);

    delete[] bits_arr;

    clear_fmpz_array(f_shares2, num_shares * num_values);

    std::cout << "Share convert time: " << sec_from(start) << std::endl;

    return sent_bytes;
}

// Batch of N (snips + num_input wire/share) validations
// Due to the nature of the final swap, both servers get the same valid array
// TODO: round split
const bool* const validate_snips(const size_t N,
                                 const size_t num_inputs,
                                 const int serverfd,
                                 const int server_num,
                                 Circuit* const * const circuit,
                                 const ClientPacket* const * const packet,
                                 const fmpz_t* const shares_p
                                 ) {
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
        // recv_unvalidated(clientfd, 1, tag);
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
            const std::string tag = get_tag(serverfd);

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
    const uint64_t max_val = 1ULL << msg.num_bits;
    const unsigned int total_inputs = msg.num_of_inputs;
    const size_t nbits[1] = {msg.num_bits};

    int cli_bytes = 0;
    for (unsigned int i = 0; i < total_inputs; i++) {
        cli_bytes += recv_in(clientfd, &share, sizeof(IntShare));
        const std::string tag(share.tag, share.tag + TAG_LENGTH);

        if (share_map.find(tag) != share_map.end()
            or share.val >= max_val)
            continue;
        share_map[tag] = share.val;

        // std::cout << "share[" << i << "] = " << share.val << std::endl;

        cli_bytes += recv_unvalidated(clientfd, tag, msg.num_bits);
    }

    std::cout << "Received " << total_inputs << " total shares" << std::endl;
    std::cout << "bytes from client: " << cli_bytes << std::endl;
    std::cout << "receive time: " << sec_from(start) << std::endl;

    if (STORE_TYPE != ot_store)
        ((CorrelatedStore*) correlated_store)->check_DaBits(total_inputs * msg.num_bits);

    start = clock_start();
    auto start2 = clock_start();

    int sent_bytes = 0;

    if (server_num == 1) {
        const size_t num_inputs = share_map.size();
        sent_bytes += send_size(serverfd, num_inputs);
        uint64_t* const shares = new uint64_t[num_inputs];
        int i = 0;
        for (const auto& share : share_map) {
            sent_bytes += send_tag(serverfd, share.first);
            shares[i] = share.second;
            i++;

            process_unvalidated(&share.first[0], msg.num_bits);
        }
        std::cout << "tag time: " << sec_from(start2) << std::endl;
        if (STORE_TYPE == validate_store) ((CorrelatedStore*) correlated_store)->check_DaBits(total_inputs * msg.num_bits);
        start2 = clock_start();
        fmpz_t* shares_p; new_fmpz_array(&shares_p, num_inputs);
        sent_bytes += share_convert(num_inputs, 1, nbits, shares, shares_p);
        std::cout << "convert time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        bool* const valid = new bool[num_inputs];
        recv_bool_batch(serverfd, valid, num_inputs);

        fmpz_t* b; new_fmpz_array(&b, 1);
        accumulate(num_inputs, 1, shares_p, valid, b);
        clear_fmpz_array(shares_p, num_inputs);
        delete[] shares;
        delete[] valid;

        std::cout << "accumulate time: " << sec_from(start2) << std::endl;

        send_fmpz(serverfd, b[0]);
        clear_fmpz_array(b, 1);
        std::cout << "total compute time: " << sec_from(start) << std::endl;
        std::cout << "sent server bytes: " << sent_bytes << std::endl;
        return RET_NO_ANS;
    } else {
        size_t num_inputs;
        recv_size(serverfd, num_inputs);
        uint64_t* const shares = new uint64_t[num_inputs];
        bool* const valid = new bool[num_inputs];

        for (unsigned int i = 0; i < num_inputs; i++) {
            const std::string tag = get_tag(serverfd);

            bool is_valid = (share_map.find(tag) != share_map.end());
            valid[i] = is_valid;
            if (!is_valid)
                continue;
            shares[i] = share_map[tag];

            process_unvalidated(tag, msg.num_bits);
        }
        std::cout << "tag time: " << sec_from(start2) << std::endl;
        if (STORE_TYPE == validate_store) ((CorrelatedStore*) correlated_store)->check_DaBits(total_inputs * msg.num_bits);
        start2 = clock_start();
        fmpz_t* shares_p; new_fmpz_array(&shares_p, num_inputs);
        sent_bytes += share_convert(num_inputs, 1, nbits, shares, shares_p);
        std::cout << "convert time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        sent_bytes += send_bool_batch(serverfd, valid, num_inputs);

        fmpz_t* a; new_fmpz_array(&a, 1);
        size_t num_valid = accumulate(num_inputs, 1, shares_p, valid, a);
        clear_fmpz_array(shares_p, num_inputs);
        delete[] shares;
        delete[] valid;

        fmpz_t b; fmpz_init(b);
        recv_fmpz(serverfd, b);
        std::cout << "accumulate time: " << sec_from(start2) << std::endl;
        std::cout << "Final valid count: " << num_valid << " / " << total_inputs << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;
        std::cout << "sent server bytes: " << sent_bytes << std::endl;
        if (num_valid < total_inputs * (1 - INVALID_THRESHOLD)) {
            std::cout << "Failing, This is less than the invalid threshold of " << INVALID_THRESHOLD << std::endl;
            return RET_INVALID;
        }

        fmpz_mod_add(b, b, a[0], mod_ctx);
        clear_fmpz_array(a, 1);
        ans = fmpz_get_ui(b);
        fmpz_clear(b);
        return RET_ANS;
    }
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
        std::string* const tag_list = new std::string[num_inputs];
        size_t idx = 0;
        for (const auto& share : share_map) {
            sent_bytes += send_tag(serverfd, share.first);
            tag_list[idx] = share.first;
            idx++;
        }
        std::cout << "tag time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        bool* const other_valid = new bool[num_inputs];
        recv_bool_batch(serverfd, other_valid, num_inputs);
        for (unsigned int i = 0; i < num_inputs; i++) {
            if (!other_valid[i])
                continue;
            b ^= share_map[tag_list[i]];
        }
        delete[] other_valid;

        send_uint64(serverfd, b);
        delete[] tag_list;
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
            const std::string tag = get_tag(serverfd);
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
        std::string* const tag_list = new std::string[num_inputs];
        size_t idx = 0;
        for (const auto& share : share_map) {
            sent_bytes += send_tag(serverfd, share.first);
            tag_list[idx] = share.first;
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
                b[j] ^= share_map[tag_list[i]][j];
        }
        delete[] other_valid;
        delete[] shares;
        delete[] tag_list;
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
            const std::string tag = get_tag(serverfd);
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
    const uint64_t max_val = 1ULL << msg.num_bits;
    const unsigned int total_inputs = msg.num_of_inputs;
    const size_t nbits[2] = {msg.num_bits, msg.num_bits * 2};
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

        // std::cout << "share[" << i << "] = " << share.val << ", " << share.val_squared << std::endl;

        if ((share_map.find(tag) != share_map.end())
            or (share.val >= max_val)
            or (share.val_squared >= max_val * max_val)
            or (packet_bytes <= 0)
            ) {
            delete packet;
            continue;
        }
        share_map[tag] = {share.val, share.val_squared, packet};

        cli_bytes += recv_unvalidated(clientfd, tag, msg.num_bits * num_dabits);
    }

    std::cout << "Received " << total_inputs << " total shares" << std::endl;
    std::cout << "bytes from client: " << cli_bytes << std::endl;
    std::cout << "receive time: " << sec_from(start) << std::endl;

    if (STORE_TYPE != ot_store)
        ((CorrelatedStore*) correlated_store)->check_DaBits(total_inputs * msg.num_bits * num_dabits);

    start = clock_start();

    auto start2 = clock_start();

    int sent_bytes = 0;

    if (server_num == 1) {
        const size_t num_inputs = share_map.size();
        sent_bytes += send_size(serverfd, num_inputs);

        uint64_t* const shares = new uint64_t[2 * num_inputs];
        ClientPacket** const packet = new ClientPacket*[num_inputs];

        std::string* const tag_list = new std::string[num_inputs];
        Circuit** const circuit = new Circuit*[num_inputs];

        size_t idx = 0;
        for (const auto& share : share_map) {
            sent_bytes += send_tag(serverfd, share.first);
            tag_list[idx] = share.first;
            idx++;

            process_unvalidated(&share.first[0], msg.num_bits * num_dabits);
        }
        std::cout << "tag time: " << sec_from(start2) << std::endl;
        if (STORE_TYPE == validate_store) ((CorrelatedStore*) correlated_store)->check_DaBits(total_inputs * msg.num_bits * num_dabits);
        start2 = clock_start();

        for (unsigned int i = 0; i < num_inputs; i++) {
            uint64_t val = 0, val2 = 0;
            std::tie(val, val2, packet[i]) = share_map[tag_list[i]];
            circuit[i] = CheckVar();
            shares[2 * i] = val;
            shares[2 * i + 1] = val2;
        }
        fmpz_t* shares_p; new_fmpz_array(&shares_p, num_inputs * 2);
        sent_bytes += share_convert(num_inputs, 2, nbits, shares, shares_p);
        std::cout << "convert time: " << sec_from(start2) << std::endl;
        start2 = clock_start();
        const bool* const snip_valid = validate_snips(
            num_inputs, 2, serverfd, server_num, circuit, packet, shares_p);

        bool* const valid = new bool[num_inputs];
        recv_bool_batch(serverfd, valid, num_inputs);

        for (unsigned int i = 0; i < num_inputs; i++) {
            valid[i] &= snip_valid[i];

            delete circuit[i];
            delete packet[i];
        }
        delete[] snip_valid;
        delete[] tag_list;
        delete[] circuit;
        delete[] packet;
        delete[] shares;
        std::cout << "validate time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        // Convert
        fmpz_t* b; new_fmpz_array(&b, 2);
        accumulate(num_inputs, 2, shares_p, valid, b);
        delete[] valid;
        clear_fmpz_array(shares_p, num_inputs * 2);

        std::cout << "accumulate time: " << sec_from(start2) << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;

        send_fmpz(serverfd, b[0]);
        send_fmpz(serverfd, b[1]);

        clear_fmpz_array(b, 2);

        std::cout << "sent non-snip server bytes: " << sent_bytes << std::endl;
        return RET_NO_ANS;
    } else {
        size_t num_inputs;
        recv_size(serverfd, num_inputs);

        uint64_t* const shares = new uint64_t[2 * num_inputs];
        ClientPacket** const packet = new ClientPacket*[num_inputs];

        bool* const valid = new bool[num_inputs];
        std::string* const tag_list = new std::string[num_inputs];
        Circuit** const circuit = new Circuit*[num_inputs];

        for (unsigned int i = 0; i < num_inputs; i++) {
            const std::string tag = get_tag(serverfd);
            tag_list[i] = tag;
            valid[i] = (share_map.find(tag) != share_map.end());

            process_unvalidated(tag, msg.num_bits * num_dabits);
        }
        std::cout << "tag time: " << sec_from(start2) << std::endl;
        if (STORE_TYPE == validate_store) ((CorrelatedStore*) correlated_store)->check_DaBits(total_inputs * msg.num_bits * num_dabits);
        start2 = clock_start();
        for (unsigned int i = 0; i < num_inputs; i++) {
            uint64_t val = 0, val2 = 0;
            if (valid[i]) {
                std::tie(val, val2, packet[i]) = share_map[tag_list[i]];
            } else {
                packet[i] = new ClientPacket(NMul);  // mock empty packet
            }
            circuit[i] = CheckVar();
            shares[2 * i] = val;
            shares[2 * i + 1] = val2;
        }
        fmpz_t* shares_p; new_fmpz_array(&shares_p, num_inputs * 2);
        sent_bytes += share_convert(num_inputs, 2, nbits, shares, shares_p);
        std::cout << "convert time: " << sec_from(start2) << std::endl;
        start2 = clock_start();
        const bool* const snip_valid = validate_snips(
            num_inputs, 2, serverfd, server_num, circuit, packet, shares_p);

        for (unsigned int i = 0; i < num_inputs; i++) {
            valid[i] &= snip_valid[i];

            delete circuit[i];
            delete packet[i];
        }
        // Send valid back, to also encapsulate pre-snip valid[]
        sent_bytes += send_bool_batch(serverfd, valid, num_inputs);
        delete[] snip_valid;
        delete[] tag_list;
        delete[] circuit;
        delete[] packet;
        delete[] shares;
        std::cout << "validate time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        // Convert
        fmpz_t* a; new_fmpz_array(&a, 2);
        size_t num_valid = accumulate(num_inputs, 2, shares_p, valid, a);

        std::cout << "accumulate time: " << sec_from(start2) << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;
        start2 = clock_start();

        delete[] valid;
        clear_fmpz_array(shares_p, num_inputs * 2);

        fmpz_t b; fmpz_init(b); recv_fmpz(serverfd, b);
        fmpz_t b2; fmpz_init(b2); recv_fmpz(serverfd, b2);

        std::cout << "Final valid count: " << num_valid << " / " << total_inputs << std::endl;
        std::cout << "sent non-snip server bytes: " << sent_bytes << std::endl;
        if (num_valid < total_inputs * (1 - INVALID_THRESHOLD)) {
            std::cout << "Failing, This is less than the invalid threshold of " << INVALID_THRESHOLD << std::endl;
            clear_fmpz_array(a, 2);
            fmpz_clear(b);
            fmpz_clear(b2);
            return RET_INVALID;
        }

        fmpz_mod_add(b, b, a[0], mod_ctx);
        fmpz_mod_add(b2, b2, a[1], mod_ctx);

        const double ex = fmpz_get_d(b) / num_valid;
        const double ex2 = fmpz_get_d(b2) / num_valid;
        ans = ex2 - (ex * ex);
        if (msg.type == VAR_OP) {
            std::cout << "Ans: " << ex2 << " - (" << ex << ")^2 = " << ans << std::endl;
        }
        if (msg.type == STDDEV_OP) {
            ans = sqrt(ans);
            std::cout << "Ans: sqrt(" << ex2 << " - (" << ex << ")^2) = " << ans << std::endl;
        }
        clear_fmpz_array(a, 2);
        fmpz_clear(b);
        fmpz_clear(b2);
        std::cout << "eval time: " << sec_from(start2) << std::endl;
        return RET_ANS;
    }
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

    const uint64_t max_val = 1ULL << msg.num_bits;
    const unsigned int total_inputs = msg.num_of_inputs;
    size_t nbits[num_fields];
    for (unsigned int i = 0; i < num_fields; i++)
        nbits[i] = msg.num_bits * (i >= degree ? 2 : 1);

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

        cli_bytes += recv_unvalidated(clientfd, tag, msg.num_bits * num_dabits);
    }

    std::cout << "Received " << total_inputs << " total shares" << std::endl;
    std::cout << "bytes from client: " << cli_bytes << std::endl;
    std::cout << "receive time: " << sec_from(start) << std::endl;

    if (STORE_TYPE != ot_store)
        ((CorrelatedStore*) correlated_store)->check_DaBits(total_inputs * msg.num_bits * num_dabits);

    start = clock_start();
    auto start2 = clock_start();

    int sent_bytes = 0;

    if (server_num == 1) {
        const size_t num_inputs = share_map.size();
        sent_bytes += send_size(serverfd, num_inputs);

        uint64_t* const shares = new uint64_t[num_inputs * num_fields];
        ClientPacket** const packet = new ClientPacket*[num_inputs];

        std::string* const tag_list = new std::string[num_inputs];
        Circuit** const circuit = new Circuit*[num_inputs];

        size_t idx = 0;
        for (const auto& share : share_map) {
            sent_bytes += send_tag(serverfd, share.first);
            tag_list[idx] = share.first;
            idx++;

            process_unvalidated(&share.first[0], msg.num_bits * num_dabits);
        }
        std::cout << "tag time: " << sec_from(start2) << std::endl;
        if (STORE_TYPE == validate_store) ((CorrelatedStore*) correlated_store)->check_DaBits(total_inputs * msg.num_bits * num_dabits);
        start2 = clock_start();

        for (unsigned int i = 0; i < num_inputs; i++) {
            uint64_t* x_vals;
            uint64_t y_val = 0;
            uint64_t* x2_vals;
            uint64_t* xy_vals;

            std::tie(x_vals, y_val, x2_vals, xy_vals, packet[i]) = share_map[tag_list[i]];
            circuit[i] = CheckLinReg(degree);

            memcpy(&shares[num_fields * i],
                   x_vals, num_x * sizeof(uint64_t));
            shares[num_fields * i + num_x] = y_val;
            memcpy(&shares[num_fields * i + num_x + 1],
                   x2_vals, num_quad * sizeof(uint64_t));
            memcpy(&shares[num_fields * i + num_x + num_quad + 1],
                   xy_vals, num_x * sizeof(uint64_t));

            delete x_vals;
            delete x2_vals;
            delete xy_vals;
        }
        fmpz_t* shares_p; new_fmpz_array(&shares_p, num_inputs * num_fields);
        sent_bytes += share_convert(num_inputs, num_fields, nbits, shares, shares_p);
        std::cout << "convert time: " << sec_from(start2) << std::endl;
        start2 = clock_start();
        const bool* const snip_valid = validate_snips(
            num_inputs, num_fields, serverfd, server_num, circuit,
            packet, shares_p);

        bool* const valid = new bool[num_inputs];
        recv_bool_batch(serverfd, valid, num_inputs);

        for (unsigned int i = 0; i < num_inputs; i++) {
            valid[i] &= snip_valid[i];

            delete circuit[i];
            delete packet[i];
        }
        delete[] snip_valid;
        delete[] tag_list;
        delete[] circuit;
        delete[] packet;
        delete[] shares;
        std::cout << "validate time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        // Convert
        fmpz_t* b; new_fmpz_array(&b, num_fields);
        accumulate(num_inputs, num_fields, shares_p, valid, b);
        delete[] valid;
        clear_fmpz_array(shares_p, num_inputs * num_fields);

        std::cout << "accumulate time: " << sec_from(start2) << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;

        send_fmpz_batch(serverfd, b, num_fields);

        clear_fmpz_array(b, num_fields);

        std::cout << "sent non-snip server bytes: " << sent_bytes << std::endl;

        return RET_NO_ANS;
    } else {
        size_t num_inputs;
        recv_size(serverfd, num_inputs);

        uint64_t* const shares = new uint64_t[num_inputs * num_fields];
        ClientPacket** const packet = new ClientPacket*[num_inputs];

        bool* const valid = new bool[num_inputs];
        std::string* const tag_list = new std::string[num_inputs];
        Circuit** const circuit = new Circuit*[num_inputs];

        for (unsigned int i = 0; i < num_inputs; i++) {
            const std::string tag = get_tag(serverfd);
            tag_list[i] = tag;
            valid[i] = (share_map.find(tag) != share_map.end());

            process_unvalidated(tag, msg.num_bits * num_dabits);
        }
        std::cout << "tag time: " << sec_from(start2) << std::endl;
        if (STORE_TYPE == validate_store) ((CorrelatedStore*) correlated_store)->check_DaBits(total_inputs * msg.num_bits * num_dabits);
        start2 = clock_start();

        for (unsigned int i = 0; i < num_inputs; i++) {
            uint64_t* x_vals;
            uint64_t y_val = 0;
            uint64_t* x2_vals;
            uint64_t* xy_vals;
            if (valid[i]) {
                std::tie(x_vals, y_val, x2_vals, xy_vals, packet[i]) = share_map[tag_list[i]];
            } else {
                x_vals = new uint64_t[num_x];
                x2_vals = new uint64_t[num_quad];
                xy_vals = new uint64_t[num_x];
                packet[i] = new ClientPacket(NMul);
            }
            circuit[i] = CheckLinReg(degree);
            memcpy(&shares[num_fields * i],
                   x_vals, num_x * sizeof(uint64_t));
            shares[num_fields * i + num_x] = y_val;
            memcpy(&shares[num_fields * i + num_x + 1],
                   x2_vals, num_quad * sizeof(uint64_t));
            memcpy(&shares[num_fields * i + num_x + num_quad + 1],
                   xy_vals, num_x * sizeof(uint64_t));

            delete x_vals;
            delete x2_vals;
            delete xy_vals;
        }
        fmpz_t* shares_p; new_fmpz_array(&shares_p, num_inputs * num_fields);
        sent_bytes += share_convert(num_inputs, num_fields, nbits, shares, shares_p);
        std::cout << "convert time: " << sec_from(start2) << std::endl;
        start2 = clock_start();
        const bool* const snip_valid = validate_snips(
            num_inputs, num_fields, serverfd, server_num, circuit,
            packet, shares_p);

        for (unsigned int i = 0; i < num_inputs; i++) {
            valid[i] &= snip_valid[i];

            delete circuit[i];
            delete packet[i];
        }
        // Send valid back, to also encapsulate pre-snip valid[]
        sent_bytes += send_bool_batch(serverfd, valid, num_inputs);
        delete[] snip_valid;
        delete[] tag_list;
        delete[] circuit;
        delete[] packet;
        delete[] shares;
        std::cout << "validate time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        // Convert
        fmpz_t* a; new_fmpz_array(&a, num_fields);
        size_t num_valid = accumulate(num_inputs, num_fields, shares_p, valid, a);

        delete[] valid;
        clear_fmpz_array(shares_p, num_inputs * num_fields);

        std::cout << "accumulate time: " << sec_from(start2) << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;
        start2 = clock_start();

        uint64_t* const x_accum = new uint64_t[degree + num_quad];
        memset(x_accum, 0, (degree + num_quad) * sizeof(uint64_t));
        uint64_t* const y_accum = new uint64_t[degree];
        memset(y_accum, 0, (degree) * sizeof(uint64_t));

        fmpz_t b; fmpz_init(b);

        x_accum[0] = num_valid;
        for (unsigned int j = 0; j < num_x; j++) {
            recv_fmpz(serverfd, b);
            fmpz_mod_add(b, b, a[j], mod_ctx);
            x_accum[1 + j] = fmpz_get_si(b);
        }
        recv_fmpz(serverfd, b);
        fmpz_mod_add(b, b, a[num_x], mod_ctx);
        y_accum[0] = fmpz_get_si(b);

        for (unsigned int j = 0; j < num_quad; j++) {
            recv_fmpz(serverfd, b);
            fmpz_mod_add(b, b, a[degree + j], mod_ctx);
            x_accum[1 + num_x + j] = fmpz_get_si(b);
        }
        for (unsigned int j = 0; j < num_x; j++) {
            recv_fmpz(serverfd, b);
            fmpz_mod_add(b, b, a[degree + num_quad + j], mod_ctx);
            y_accum[1 + j] = fmpz_get_si(b);
        }

        clear_fmpz_array(a, num_fields);
        fmpz_clear(b);

        std::cout << "Final valid count: " << num_valid << " / " << total_inputs << std::endl;
        std::cout << "sent non-snip server bytes: " << sent_bytes << std::endl;
        if (num_valid < total_inputs * (1 - INVALID_THRESHOLD)) {
            std::cout << "Failing, This is less than the invalid threshold of " << INVALID_THRESHOLD << std::endl;
            delete[] x_accum;
            delete[] y_accum;
            return RET_INVALID;
        }

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
        start2 = clock_start();

        bool* const valid = new bool[num_inputs];
        fmpz_t* sums; new_fmpz_array(&sums, num_inputs);
        bool total_parity = 0;
        fmpz_t sum; fmpz_init_set_si(sum, 0);
        for (unsigned int i = 0; i < num_inputs; i++) {
            fmpz_zero(sums[i]);
            for (unsigned int j = 0; j < max_inp; j++) {
                // std::cout << "share_p[" << i << ", " << j << "] = ";
                // fmpz_print(shares_p[i * max_inp + j]); std::cout << std::endl;

                fmpz_mod_add(sums[i], sums[i], shares_p[i * max_inp + j], mod_ctx);
                total_parity ^= shares[i * max_inp + j];
            }
            fmpz_mod_add(sum, sum, sums[i], mod_ctx);
        }
        // Batch check.
        bool total_parity_other;
        fmpz_t sum_other; fmpz_init(sum_other);
        sent_bytes += send_bool(serverfd, total_parity);
        sent_bytes += send_fmpz(serverfd, sum);
        recv_bool(serverfd, total_parity_other);
        recv_fmpz(serverfd, sum_other);
        fmpz_clear(sum);
        delete[] shares;

        bool all_valid = (total_parity ^ total_parity_other) == (total_inputs % 2);
        if (all_valid) {
            fmpz_mod_add(sum, sum, sum_other, mod_ctx);
            all_valid = fmpz_equal_ui(sum, total_inputs);
        }
        if (all_valid) {
            memset(valid, true, num_inputs);
        } else {
            // TODO: Single check (or binary search)
            memset(valid, false, num_inputs);
            std::cout << "Batch not valid. Individual check currently not implemented" << std::endl;
        }

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
            const std::string tag = get_tag(serverfd);
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
        start2 = clock_start();

        // Validity check
        bool* const parity = new bool[num_inputs];
        bool total_parity = 0;
        fmpz_t* sums; new_fmpz_array(&sums, num_inputs);
        fmpz_t sum; fmpz_init_set_si(sum, 0);
        for (unsigned int i = 0; i < num_inputs; i++) {
            parity[i] = false;
            fmpz_zero(sums[i]);
            for (unsigned int j = 0; j < max_inp; j++) {
                // std::cout << "share_p[" << i << ", " << j << "] = ";
                // fmpz_print(shares_p[i * max_inp + j]); std::cout << std::endl;

                parity[i] ^= shares[i * max_inp + j];
                fmpz_mod_add(sums[i], sums[i], shares_p[i * max_inp + j], mod_ctx);
            }
            total_parity ^= parity[i];
            fmpz_mod_add(sum, sum, sums[i], mod_ctx);
        }
        // batch check
        bool total_parity_other;
        fmpz_t sum_other; fmpz_init(sum_other);
        sent_bytes += send_bool(serverfd, total_parity);
        sent_bytes += send_fmpz(serverfd, sum);
        recv_bool(serverfd, total_parity_other);
        recv_fmpz(serverfd, sum_other);
        bool all_valid = false;
        all_valid = (total_parity ^ total_parity_other) == (total_inputs % 2);
        if (all_valid) {
            fmpz_mod_add(sum, sum, sum_other, mod_ctx);
            all_valid = fmpz_equal_ui(sum, total_inputs);
        }
        sent_bytes += send_bool(serverfd, all_valid);
        if (all_valid) {
            memset(valid, true, num_inputs);
        } else {
            // Batch check fails. Test single. Binary search is more rounds but (possibly) less bytes).

            bool* const parity_other = new bool[num_inputs];
            fmpz_t* sums_other; new_fmpz_array(&sums_other, num_inputs);
            recv_bool_batch(serverfd, parity_other, num_inputs);
            recv_fmpz_batch(serverfd, sums_other, num_inputs);
            for (unsigned int i = 0; i < num_inputs; i++) {
                if (!valid[i])
                    continue;
                if ((parity[i] ^ parity_other[i]) == 0) {
                    valid[i] = false;
                    continue;
                }
                fmpz_mod_add(sums[i], sums[i], sums_other[i], mod_ctx);
                if (!fmpz_is_one(sums[i]))
                    valid[i] = false;
            }
            clear_fmpz_array(sums_other, num_inputs);
            delete[] parity_other;

            sent_bytes += send_bool_batch(serverfd, valid, num_inputs);
        }

        clear_fmpz_array(sums, num_inputs);
        delete[] parity;
        std::cout << "validate time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        fmpz_t* a; new_fmpz_array(&a, max_inp);
        size_t num_valid = accumulate(num_inputs, max_inp, shares_p, valid, a);
        delete[] valid;
        delete[] shares;
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

    // Not actually necessary for evaluation
    // Could help for "verify" that hashes actually have right sign after finding the output, but it doesn't help evaluation.
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
    if (server_num == 1) {
        const size_t num_inputs = share_map.size();
        const size_t n = num_inputs * b;
        sent_bytes += send_size(serverfd, num_inputs);
        bool* const valid = new bool[num_inputs];
        bool* const x = new bool[n];
        bool* const y = new bool[n];

        size_t idx = 0;
        for (const auto& share : share_map) {
            sent_bytes += send_tag(serverfd, share.first);

            memcpy(&x[idx*b], share.second, b);
            memcpy(&y[idx*b], &(share.second[b]), b);

            // std::cout << "share[" << idx << "] = ";
            // for (unsigned int j = 0; j < b; j++)
            //     std::cout << x[idx * b + j];
            // std::cout << ", ";
            // for (unsigned int j = 0; j < b; j++)
            //     std::cout << y[idx * b + j];
            // std::cout << std::endl;

            idx++;

            delete[] share.second;
        }
        recv_bool_batch(serverfd, valid, num_inputs);
        std::cout << "tag time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        sent_bytes += correlated_store->heavy_convert(num_inputs, b, x, y, valid, bucket0, bucket1);
        delete[] x;
        delete[] y;
        std::cout << "convert+accum time: " << sec_from(start2) << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;
        std::cout << "compute bytes sent: " << sent_bytes << std::endl;

        // Evaluate
        // True means |1| is larger than |0|
        start2 = clock_start();
        fmpz_t* larger; new_fmpz_array(&larger, b);
        sent_bytes = correlated_store->abs_cmp(b, bucket0, bucket1, larger);
        sent_bytes += send_fmpz_batch(serverfd, larger, b);
        std::cout << "evaluate time: " << sec_from(start2) << std::endl;
        std::cout << "evaluate bytes sent: " << sent_bytes << std::endl;

        clear_fmpz_array(bucket0, b);
        clear_fmpz_array(bucket1, b);
        // Other side has it, no more eval needed.

        clear_fmpz_array(larger, b);
        delete[] valid;

        return RET_NO_ANS;

    } else {
        size_t num_inputs;
        recv_size(serverfd, num_inputs);
        const size_t n = num_inputs * b;
        bool* const x = new bool[n];
        bool* const y = new bool[n];
        bool* const valid = new bool[num_inputs];

        bool* share;
        for (unsigned int i = 0; i < num_inputs; i++) {
            const std::string tag = get_tag(serverfd);
            valid[i] = (share_map.find(tag) != share_map.end());
            if (!valid[i]) {
                memset(&x[i * b], 0, b);
                memset(&y[i * b], 0, b);
                continue;
            }
            share = share_map[tag];
            memcpy(&x[i * b], share, b);
            memcpy(&y[i * b], &(share[b]), b);

            // std::cout << "share[" << i << "] = ";
            // for (unsigned int j = 0; j < b; j++)
            //     std::cout << x[i * b + j];
            // std::cout << ", ";
            // for (unsigned int j = 0; j < b; j++)
            //     std::cout << y[i * b + j];
            // std::cout << std::endl;

            delete[] share;
        }
        sent_bytes += send_bool_batch(serverfd, valid, num_inputs);
        std::cout << "tag time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        sent_bytes += correlated_store->heavy_convert(num_inputs, b, x, y, valid, bucket0, bucket1);
        delete[] x;
        delete[] y;
        std::cout << "convert+accum time: " << sec_from(start2) << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;
        std::cout << "compute bytes sent: " << sent_bytes << std::endl;

        // Evaluate

        start2 = clock_start();
        fmpz_t* larger_0; new_fmpz_array(&larger_0, b);
        fmpz_t* larger_1; new_fmpz_array(&larger_1, b);
        sent_bytes = correlated_store->abs_cmp(b, bucket0, bucket1, larger_0);
        recv_fmpz_batch(serverfd, larger_1, b);

        clear_fmpz_array(bucket0, b);
        clear_fmpz_array(bucket1, b);

        size_t num_valid = 0;
        for (unsigned int i = 0; i < num_inputs; i++)
            num_valid += valid[i];
        std::cout << "Final valid count: " << num_valid << " / " << total_inputs << std::endl;
        if (num_valid < total_inputs * (1 - INVALID_THRESHOLD)) {
            std::cout << "Failing, This is less than the invalid threshold of " << INVALID_THRESHOLD << std::endl;
            clear_fmpz_array(larger_0, b);
            clear_fmpz_array(larger_1, b);
            delete[] valid;
            return RET_INVALID;
        }

        // fmpz_from_bool_array?
        uint64_t ans = 0;
        bool* larger = new bool[b];
        for (unsigned int j = 0; j < b; j++) {
            fmpz_mod_add(larger_0[j], larger_0[j], larger_1[j], mod_ctx);
            larger[j] = fmpz_is_one(larger_0[j]);
            // std::cout << "bucket" << (larger[j] ? 1 : 0) << "[" << j << "] is heavier" << std::endl;
            ans |= (larger[j] << j);
        }

        std::cout << "evaluate time: " << sec_from(start2) << std::endl;
        std::cout << "evaluate sent bytes: " << sent_bytes << std::endl;

        std::cout << "### Heavy hitter value is " << ans << std::endl;

        clear_fmpz_array(larger_0, b);
        clear_fmpz_array(larger_1, b);
        delete[] larger;
        delete[] valid;
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
    const size_t share_size_sh = cfg.Q * cfg.SH_depth;
    // For each Q, identify b (first hash)
    const size_t share_size_mask = cfg.Q * cfg.B;
    // Count-min
    const size_t share_size_count = cfg.countmin_cfg.d * cfg.countmin_cfg.w;
    // SingleHeavy bucket pairs
    const size_t num_sh = cfg.Q * cfg.B * cfg.SH_depth;
    // Values being converted in B2A
    const size_t share_convert_size = share_size_sh + share_size_count + share_size_mask;

    if (STORE_TYPE != ot_store) {
        ((CorrelatedStore*) correlated_store)->check_DaBits(total_inputs * share_convert_size);
    }

    fmpz_t* bucket0; new_fmpz_array(&bucket0, num_sh);
    fmpz_t* bucket1; new_fmpz_array(&bucket1, num_sh);

    auto start2 = clock_start();
    for (unsigned int i = 0; i < total_inputs; i++) {
        char tag_c[TAG_LENGTH];
        cli_bytes += recv_in(clientfd, &tag_c[0], TAG_LENGTH);
        const std::string tag(tag_c, tag_c+TAG_LENGTH);

        bool* const share_sh_x = new bool[share_size_sh];
        bool* const share_sh_y = new bool[share_size_sh];
        bool* const share_mask = new bool[share_size_mask];
        bool* const share_count = new bool[share_size_count];
        cli_bytes += recv_bool_batch(clientfd, share_sh_x, share_size_sh);
        cli_bytes += recv_bool_batch(clientfd, share_sh_y, share_size_sh);
        cli_bytes += recv_bool_batch(clientfd, share_mask, share_size_mask);
        cli_bytes += recv_bool_batch(clientfd, share_count, share_size_count);

        if (share_map.find(tag) != share_map.end()) {
            delete[] share_sh_x;
            delete[] share_sh_y;
            delete[] share_mask;
            delete[] share_count;
            continue;
        }

        share_map[tag] = {share_sh_x, share_sh_y, share_mask, share_count};
    }

    std::cout << "Received " << total_inputs << " total shares" << std::endl;
    std::cout << "bytes from client: " << cli_bytes << std::endl;
    std::cout << "receive time: " << sec_from(start2) << std::endl;

    /* Stages:
    Validate each b (share_mask) is frequency vector. With B2A to check sum to 1
    "Validate H_q vector correct". share_sh. Not needed for R=1, is just entry
    (joint share validity)
    Update Countmin with share_count (and freq check in parallel)
    For each (q, b):
      final mask = share_sh & R_mask.
      SH update, except OT multiply [z] by mask.

    TODO: round collapse

    mask: B2A just for validation. Kept as bool for Ot multiply
    sh: half B2A for use (z), half kept selection
    count: B2A for validation and accumulation
    */

    start = clock_start();
    start2 = clock_start();
    int64_t sent_bytes = 0;
    if (server_num == 1) {
        const size_t num_inputs = share_map.size();
        sent_bytes += send_size(serverfd, num_inputs);
        std::cout << "num_inputs: " << num_inputs << std::endl;

        std::string* const tag_list = new std::string[num_inputs];
        bool* const valid = new bool[num_inputs];
        bool* const shares_sh_x = new bool[num_inputs * share_size_sh];
        bool* const shares_sh_y = new bool[num_inputs * share_size_sh];
        bool* const shares_mask = new bool[num_inputs * share_size_mask];
        bool* const shares_count = new bool[num_inputs * share_size_count];

        size_t idx = 0;
        for (const auto& share : share_map) {
            sent_bytes += send_tag(serverfd, share.first);
            tag_list[idx] = share.first;
            idx++;
        }
        recv_bool_batch(serverfd, valid, num_inputs);
        std::cout << "tag time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        for (unsigned int i = 0; i < num_inputs; i++) {
            if (!valid[i]) continue;
            bool* a; bool* b; bool* c; bool* d;
            std::tie(a, b, c, d) = share_map[tag_list[i]];
            memcpy(&shares_sh_x[i * share_size_sh], a, share_size_sh);
            memcpy(&shares_sh_y[i * share_size_sh], b, share_size_sh);
            memcpy(&shares_mask[i * share_size_mask], c, share_size_mask);
            memcpy(&shares_count[i * share_size_count], d, share_size_count);
            delete a;
            delete b;
            delete c;
            delete d;
        }
        delete[] tag_list;

        auto start3 = clock_start();

        // Single round all B2A: count, mask, y
        const size_t convert_size = num_inputs * share_convert_size;
        fmpz_t* shares_p; new_fmpz_array(&shares_p, convert_size);
        bool* shares_2 = new bool[convert_size];
        memcpy(shares_2, shares_count, num_inputs * share_size_count);
        delete[] shares_count;
        memcpy(&shares_2[num_inputs * share_size_count],
                shares_mask, num_inputs * share_size_mask);
        const size_t share_y_offset = num_inputs * (share_size_count + share_size_mask);
        memcpy(&shares_2[share_y_offset], shares_sh_y, num_inputs * share_size_sh);
        delete[] shares_sh_y;
        sent_bytes += correlated_store->b2a_single(convert_size, shares_2, shares_p);
        // We use first part of shares_p as the countmin shares.

        // Freq check. Batch across all inputs
        // count: num * d entries size w
        // mask: num * Q entries size mask
        bool* parity = new bool[cfg.Q + cfg.countmin_cfg.d];
        fmpz_t* sums; new_fmpz_array(&sums, cfg.Q + cfg.countmin_cfg.d);
        idx = 0;
        for (unsigned int i = 0; i < num_inputs; i++) {
            for (unsigned int j = 0; j < cfg.countmin_cfg.d; j++) {
                for (unsigned int k = 0; k < cfg.countmin_cfg.w; k++) {
                    parity[j + cfg.Q] ^= shares_2[idx];
                    fmpz_mod_add(sums[j], sums[j], shares_p[idx], mod_ctx);
                    idx++;
                }
            }
        }
        for (unsigned int i = 0; i < num_inputs; i++) {
            for (unsigned int j = 0; j < cfg.Q; j++) {
                for (unsigned int k = 0; k < cfg.B; k++) {
                    parity[j] ^= shares_2[idx];
                    fmpz_mod_add(sums[j + cfg.countmin_cfg.d], sums[j+ cfg.countmin_cfg.d], shares_p[idx], mod_ctx);
                    idx++;
                }
            }
        }
        delete[] shares_2;
        bool* parity_other = new bool[cfg.Q + cfg.countmin_cfg.d];
        fmpz_t* sums_other; new_fmpz_array(&sums_other, cfg.Q + cfg.countmin_cfg.d);
        sent_bytes += send_bool_batch(serverfd, parity, cfg.Q + cfg.countmin_cfg.d);
        sent_bytes += send_fmpz_batch(serverfd, sums, cfg.Q + cfg.countmin_cfg.d);
        recv_bool_batch(serverfd, parity_other, cfg.Q + cfg.countmin_cfg.d);
        recv_fmpz_batch(serverfd, sums_other, cfg.Q + cfg.countmin_cfg.d);
        bool all_valid = false;
        for (unsigned int j = 0; j < cfg.Q + cfg.countmin_cfg.d; j++) {
            all_valid |= (parity[j] ^ parity_other[j]) == (total_inputs % 2);
            fmpz_mod_add(sums[j], sums[j], sums_other[j], mod_ctx);
            all_valid |= fmpz_equal_ui(sums[j], total_inputs);
            if (!all_valid)
                continue;
        }
        delete[] parity;
        delete[] parity_other;
        clear_fmpz_array(sums, cfg.Q + cfg.countmin_cfg.d);
        clear_fmpz_array(sums_other, cfg.Q + cfg.countmin_cfg.d);
        if (all_valid) {
            memset(valid, true, num_inputs);
        } else {
            memset(valid, false, num_inputs);
            std::cout << "Batch not valid. Individual check currently not implemented" << std::endl;
        }

        std::cout << "first B2A and freq validation time: " << sec_from(start3) << std::endl;

        fmpz_t* countmin_accum; new_fmpz_array(&countmin_accum, share_size_count);
        accumulate(num_inputs, share_size_count, shares_p, valid, countmin_accum);

        start3 = clock_start();
        sent_bytes += correlated_store->heavy_convert_mask(
            num_inputs, cfg.Q, cfg.B, cfg.SH_depth,
            shares_sh_x, &shares_p[share_y_offset], shares_mask,
            valid, bucket0, bucket1);
        clear_fmpz_array(shares_p, convert_size);
        delete[] shares_sh_x;
        delete[] shares_mask;
        std::cout << "heavy_convert time: " << sec_from(start3) << std::endl;
        std::cout << "convert+accum time: " << sec_from(start2) << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;
        std::cout << "compute bytes sent: " << sent_bytes << std::endl;

        start2 = clock_start();

        size_t num_valid = 0;
        for (unsigned int i = 0; i < num_inputs; i++)
            num_valid += valid[i];
        std::cout << "Final valid count: " << num_valid << " / " << total_inputs << std::endl;
        delete[] valid;
        if (num_valid < total_inputs * (1 - INVALID_THRESHOLD)) {
            std::cout << "Failing, This is less than the invalid threshold of " << INVALID_THRESHOLD << std::endl;
            return RET_INVALID;
        }

        uint64_t* top_values = new uint64_t[K];
        uint64_t* top_freqs = new uint64_t[K];

        full_heavy_extract(server_num, cfg, bucket0, bucket1, hash_seed_split, hash_seed_count,
                countmin_accum, num_inputs, top_values, top_freqs);
        garbleIO->flush();

        clear_fmpz_array(bucket0, num_sh);
        clear_fmpz_array(bucket1, num_sh);

        std::cout << "Top K = " << K << " values and freqs, decreasing\n";
        for (unsigned int i = 0; i < K; i++) {
          std::cout << "Value: " << top_values[i] << ", freq: " << top_freqs[i] << std::endl;
        }

        delete[] top_values;
        delete[] top_freqs;

        return RET_ANS;
    } else {
        size_t num_inputs;
        recv_size(serverfd, num_inputs);
        std::cout << "num_inputs: " << num_inputs << std::endl;

        std::string* const tag_list = new std::string[num_inputs];
        bool* const valid = new bool[num_inputs];
        bool* const shares_sh_x = new bool[num_inputs * share_size_sh];
        bool* const shares_sh_y = new bool[num_inputs * share_size_sh];
        bool* const shares_mask = new bool[num_inputs * share_size_mask];
        bool* const shares_count = new bool[num_inputs * share_size_count];

        for (unsigned int i = 0; i < num_inputs; i++) {
            const std::string tag = get_tag(serverfd);
            tag_list[i] = tag;
            valid[i] = (share_map.find(tag) != share_map.end());
        }
        send_bool_batch(serverfd, valid, num_inputs);
        std::cout << "tag time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        for (unsigned int i = 0; i < num_inputs; i++) {
            if (!valid[i]) continue;
            bool* a; bool* b; bool* c; bool* d;
            std::tie(a, b, c, d) = share_map[tag_list[i]];
            memcpy(&shares_sh_x[i * share_size_sh], a, share_size_sh);
            memcpy(&shares_sh_y[i * share_size_sh], b, share_size_sh);
            memcpy(&shares_mask[i * share_size_mask], c, share_size_mask);
            memcpy(&shares_count[i * share_size_count], d, share_size_count);
            delete a;
            delete b;
            delete c;
            delete d;
        }
        delete[] tag_list;

        auto start3 = clock_start();

        // Single round all B2A: count, mask, y
        const size_t convert_size = num_inputs * share_convert_size;
        fmpz_t* shares_p; new_fmpz_array(&shares_p, convert_size);
        bool* shares_2 = new bool[convert_size];
        memcpy(shares_2, shares_count, num_inputs * share_size_count);
        memcpy(&shares_2[num_inputs * share_size_count],
                shares_mask, num_inputs * share_size_mask);
        const size_t share_y_offset = num_inputs * (share_size_count + share_size_mask);
        memcpy(&shares_2[share_y_offset], shares_sh_y, num_inputs * share_size_sh);
        delete[] shares_sh_y;
        delete[] shares_count;
        sent_bytes += correlated_store->b2a_single(convert_size, shares_2, shares_p);
        // We just use first part of shares_p as the countmin shares.

        // Freq check. Batch across all inputs
        // mask: num * Q entries size mask
        // count: num * d entries size w
        bool* parity = new bool[cfg.Q + cfg.countmin_cfg.d];
        fmpz_t* sums; new_fmpz_array(&sums, cfg.Q + cfg.countmin_cfg.d);
        size_t idx = 0;
        for (unsigned int i = 0; i < num_inputs; i++) {
            for (unsigned int j = 0; j < cfg.countmin_cfg.d; j++) {
                for (unsigned int k = 0; k < cfg.countmin_cfg.w; k++) {
                    parity[j + cfg.Q] ^= shares_2[idx];
                    fmpz_mod_add(sums[j], sums[j], shares_p[idx], mod_ctx);
                    idx++;
                }
            }
        }
        for (unsigned int i = 0; i < num_inputs; i++) {
            for (unsigned int j = 0; j < cfg.Q; j++) {
                for (unsigned int k = 0; k < cfg.B; k++) {
                    parity[j] ^= shares_2[idx];
                    fmpz_mod_add(sums[j + cfg.countmin_cfg.d], sums[j+ cfg.countmin_cfg.d], shares_p[idx], mod_ctx);
                    idx++;
                }
            }
        }
        delete[] shares_2;
        bool* parity_other = new bool[cfg.Q + cfg.countmin_cfg.d];
        fmpz_t* sums_other; new_fmpz_array(&sums_other, cfg.Q + cfg.countmin_cfg.d);
        sent_bytes += send_bool_batch(serverfd, parity, cfg.Q + cfg.countmin_cfg.d);
        sent_bytes += send_fmpz_batch(serverfd, sums, cfg.Q + cfg.countmin_cfg.d);
        recv_bool_batch(serverfd, parity_other, cfg.Q + cfg.countmin_cfg.d);
        recv_fmpz_batch(serverfd, sums_other, cfg.Q + cfg.countmin_cfg.d);
        bool all_valid = false;
        for (unsigned int j = 0; j < cfg.Q + cfg.countmin_cfg.d; j++) {
            all_valid |= (parity[j] ^ parity_other[j]) == (total_inputs % 2);
            fmpz_mod_add(sums[j], sums[j], sums_other[j], mod_ctx);
            all_valid |= fmpz_equal_ui(sums[j], total_inputs);
            if (!all_valid)
                continue;
        }
        delete[] parity;
        delete[] parity_other;
        clear_fmpz_array(sums, cfg.Q + cfg.countmin_cfg.d);
        clear_fmpz_array(sums_other, cfg.Q + cfg.countmin_cfg.d);
        if (all_valid) {
            memset(valid, true, num_inputs);
        } else {
            memset(valid, false, num_inputs);
            std::cout << "Batch not valid. Individual check currently not implemented" << std::endl;
        }

        std::cout << "first B2A and freq validation time: " << sec_from(start3) << std::endl;

        fmpz_t* countmin_accum; new_fmpz_array(&countmin_accum, share_size_count);
        accumulate(num_inputs, share_size_count, shares_p, valid, countmin_accum);

        start3 = clock_start();
        sent_bytes += correlated_store->heavy_convert_mask(
            num_inputs, cfg.Q, cfg.B, cfg.SH_depth,
            shares_sh_x, &shares_p[share_y_offset], shares_mask,
            valid, bucket0, bucket1);
        clear_fmpz_array(shares_p, convert_size);
        delete[] shares_sh_x;
        delete[] shares_mask;
        std::cout << "heavy_convert time: " << sec_from(start3) << std::endl;
        std::cout << "convert+accum time: " << sec_from(start2) << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;
        std::cout << "compute bytes sent: " << sent_bytes << std::endl;

        start2 = clock_start();

        size_t num_valid = 0;
        for (unsigned int i = 0; i < num_inputs; i++)
            num_valid += valid[i];
        std::cout << "Final valid count: " << num_valid << " / " << total_inputs << std::endl;
        delete[] valid;
        if (num_valid < total_inputs * (1 - INVALID_THRESHOLD)) {
            std::cout << "Failing, This is less than the invalid threshold of " << INVALID_THRESHOLD << std::endl;
            return RET_INVALID;
        }

        uint64_t* top_values = new uint64_t[K];
        uint64_t* top_freqs = new uint64_t[K];

        full_heavy_extract(server_num, cfg, bucket0, bucket1, hash_seed_split, hash_seed_count,
                countmin_accum, num_inputs, top_values, top_freqs);
        garbleIO->flush();

        // Extract clobbers accum?

        // reveal_fmpz_batch(serverfd, countmin_accum, share_size_count);

        clear_fmpz_array(bucket0, num_sh);
        clear_fmpz_array(bucket1, num_sh);
        // Count-min cleared elsewhere?

        std::cout << "Top K = " << K << " values and freqs, decreasing\n";
        for (unsigned int i = 0; i < K; i++) {
          std::cout << "Value: " << top_values[i] << ", freq: " << top_freqs[i] << std::endl;
        }

        delete[] top_values;
        delete[] top_freqs;

        std::cout << "eval time: " << sec_from(start2) << std::endl;

        return RET_ANS;
    }
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
    const size_t share_size_sh = cfg.Q * cfg.SH_depth;
    // For each Q, identify b (first hash)
    const size_t share_size_bucket = cfg.Q * cfg.B;
    // Layer selector.
    const size_t share_size_layer = cfg.R;
    // Count-min
    const size_t share_size_count = cfg.countmin_cfg.d * cfg.countmin_cfg.w;
    // SingleHeavy bucket pairs
    const size_t num_sh = cfg.R * cfg.Q * cfg.B * cfg.SH_depth;
    // Values being converted in B2A: count, bucket, y
    const size_t share_convert_size = share_size_count + share_size_bucket + share_size_sh;

    if (STORE_TYPE != ot_store) {
        ((CorrelatedStore*) correlated_store)->check_DaBits(total_inputs * share_convert_size);
        ((CorrelatedStore*) correlated_store)->check_BoolTriples(2 * total_inputs * share_size_bucket * share_size_layer);
    }

    fmpz_t* bucket0; new_fmpz_array(&bucket0, num_sh);
    fmpz_t* bucket1; new_fmpz_array(&bucket1, num_sh);

    // Read
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

    /* Cross product bucket x layer first:
        boolean AND
        then have that as mask into buckets (bool * arith)
        since OT for bool * arith, vs one triple for AND
        can do in parallel with B2A
    */

    /*
    Stages:
    - Validate bucket (freq), sum to 1 via B2A
    - TODO: validate layer?
    - update count-min with share_count
    - cross product AND: bucket[qb] and layer[r]
    - for each (r, q, b)
        - mask = bucket[qb] & layer[r]
        - SH update, with OT mult [z] by mask
    Rounds:
    (TODO: fit in b2a rerun, if share validation fails)
    TODO: round collapse
    TODO: shuffle timing printing
    1: B2A y (use), bucket (validation), count-min (agg, valid)
        cross multiply bucket and layer (AND)
    2: frequency checks: bucket and count-min
        Snips?
    2.5:  Heavy convert:
        SH update with [z] * mask?
        OT stuff, so doesn't overlap easily
        Also since lots of OT, by far the slowest part
        double check alignment
    Final: Accumulation (after all validation)
    */
    start = clock_start();
    start2 = clock_start();
    int64_t sent_bytes = 0;
    if (server_num == 1) {
        const size_t num_inputs = share_map.size();
        sent_bytes += send_size(serverfd, num_inputs);
        std::cout << "num_inputs: " << num_inputs << std::endl;

        std::string* const tag_list = new std::string[num_inputs];
        bool* const valid = new bool[num_inputs];
        bool* const shares_sh_x = new bool[num_inputs * share_size_sh];
        bool* const shares_sh_y = new bool[num_inputs * share_size_sh];
        bool* const shares_bucket = new bool[num_inputs * share_size_bucket];
        bool* const shares_layer = new bool[num_inputs * share_size_layer];
        bool* const shares_count = new bool[num_inputs * share_size_count];

        // Align by tag
        size_t idx = 0;
        for (const auto& share : share_map) {
            sent_bytes += send_tag(serverfd, share.first);
            tag_list[idx] = share.first;
            idx++;
        }
        recv_bool_batch(serverfd, valid, num_inputs);
        std::cout << "tag time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        for (unsigned int i = 0; i < num_inputs; i++) {
            if (!valid[i]) continue;
            bool* a; bool* b; bool* c; bool* d; bool* e;
            std::tie(a, b, c, d, e) = share_map[tag_list[i]];
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
        delete[] tag_list;

        auto start3 = clock_start();

        // Round 1: All B2A: (count, bucket, y)
        // count used (and freq checked), so first
        // bucket similarly just freq checked, so second
        // y final, with offset
        const size_t convert_size = num_inputs * share_convert_size;
        fmpz_t* shares_p; new_fmpz_array(&shares_p, convert_size);
        bool* shares_2 = new bool[convert_size];
        // Count first
        memcpy(shares_2, shares_count, num_inputs * share_size_count);
        delete[] shares_count;
        // Bucket second
        memcpy(&shares_2[num_inputs * share_size_count],
                shares_bucket, num_inputs * share_size_bucket);
        // Y third
        const size_t share_y_offset = num_inputs * (share_size_count + share_size_bucket);
        memcpy(&shares_2[share_y_offset], shares_sh_y, num_inputs * share_size_sh);
        delete[] shares_sh_y;
        sent_bytes += correlated_store->b2a_single(convert_size, shares_2, shares_p);
        std::cout << "B2A time: " << sec_from(start3) << std::endl; start3 = clock_start();

        // Round 1.5: Mask multiply (AND)
        // TODO: fold into above
        bool* mask = new bool[num_inputs * share_size_bucket * share_size_layer];
        sent_bytes += correlated_store->multiply_BoolShares_cross(
                num_inputs, share_size_bucket, share_size_layer,
                shares_bucket, shares_layer, mask);
        // Fold R into Q (just "more copies")
        std::cout << "cross time: " << sec_from(start3) << std::endl; start3 = clock_start();

        // Round 2: Frequency checks
        // Ensure count and bucket are freq vectors (one send)
        // odd parity, sum to 1
        bool* parity = new bool[cfg.Q + cfg.countmin_cfg.d];
        fmpz_t* sums; new_fmpz_array(&sums, cfg.Q + cfg.countmin_cfg.d);
        idx = 0;
        for (unsigned int i = 0; i < num_inputs; i++) {
            for (unsigned int j = 0; j < cfg.countmin_cfg.d; j++) {
                for (unsigned int k = 0; k < cfg.countmin_cfg.w; k++) {
                    parity[j + cfg.Q] ^= shares_2[idx];
                    fmpz_mod_add(sums[j], sums[j], shares_p[idx], mod_ctx);
                    idx++;
                }
            }
        }
        for (unsigned int i = 0; i < num_inputs; i++) {
            for (unsigned int j = 0; j < cfg.Q; j++) {
                for (unsigned int k = 0; k < cfg.B; k++) {
                    parity[j] ^= shares_2[idx];
                    fmpz_mod_add(sums[j + cfg.countmin_cfg.d], sums[j+ cfg.countmin_cfg.d], shares_p[idx], mod_ctx);
                    idx++;
                }
            }
        }
        delete[] shares_2;
        bool* parity_other = new bool[cfg.Q + cfg.countmin_cfg.d];
        fmpz_t* sums_other; new_fmpz_array(&sums_other, cfg.Q + cfg.countmin_cfg.d);

        // Not reveal wrapper, since paired together.
        sent_bytes += send_bool_batch(serverfd, parity, cfg.Q + cfg.countmin_cfg.d);
        sent_bytes += send_fmpz_batch(serverfd, sums, cfg.Q + cfg.countmin_cfg.d);
        recv_bool_batch(serverfd, parity_other, cfg.Q + cfg.countmin_cfg.d);
        recv_fmpz_batch(serverfd, sums_other, cfg.Q + cfg.countmin_cfg.d);

        bool all_valid = false;
        for (unsigned int j = 0; j < cfg.Q + cfg.countmin_cfg.d; j++) {
            all_valid |= (parity[j] ^ parity_other[j]) == (total_inputs % 2);
            fmpz_mod_add(sums[j], sums[j], sums_other[j], mod_ctx);
            all_valid |= fmpz_equal_ui(sums[j], total_inputs);
            if (!all_valid)
                continue;
        }
        delete[] parity;
        delete[] parity_other;
        clear_fmpz_array(sums, cfg.Q + cfg.countmin_cfg.d);
        clear_fmpz_array(sums_other, cfg.Q + cfg.countmin_cfg.d);
        if (all_valid) {
            memset(valid, true, num_inputs);
        } else {
            memset(valid, false, num_inputs);
            std::cout << "Batch not valid. Individual check currently not implemented" << std::endl;
        }
        std::cout << "freq validation time: " << sec_from(start3) << std::endl;


        // Round 2.5: Heavy convert with mask
        // 2 rounds of OT, which can't easily fold
        // Also accumulates into buckets
        start3 = clock_start();
        sent_bytes += correlated_store->heavy_convert_mask(
            num_inputs, cfg.Q, cfg.B * cfg.R, cfg.SH_depth,
            shares_sh_x, &shares_p[share_y_offset], mask, valid, bucket0, bucket1);
        // Also count-min accumulation
        fmpz_t* countmin_accum; new_fmpz_array(&countmin_accum, share_size_count);
        accumulate(num_inputs, share_size_count, shares_p, valid, countmin_accum);
        clear_fmpz_array(shares_p, convert_size);
        delete[] shares_sh_x;
        delete[] mask;
        std::cout << "heavy_convert time: " << sec_from(start3) << std::endl;
        std::cout << "convert+accum time: " << sec_from(start2) << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;
        std::cout << "compute bytes sent: " << sent_bytes << std::endl;

        start2 = clock_start();

        // Evaluation: Validity check
        size_t num_valid = 0;
        for (unsigned int i = 0; i < num_inputs; i++)
            num_valid += valid[i];
        std::cout << "Final valid count: " << num_valid << " / " << total_inputs << std::endl;
        delete[] valid;
        if (num_valid < total_inputs * (1 - INVALID_THRESHOLD)) {
            std::cout << "Failing, This is less than the invalid threshold of " << INVALID_THRESHOLD << std::endl;
            return RET_INVALID;
        }

        // Evaluation: Answer extraction
        start2 = clock_start();

        uint64_t* top_values = new uint64_t[K];
        uint64_t* top_freqs = new uint64_t[K];

        full_heavy_extract(server_num, cfg, bucket0, bucket1, hash_seed_split, hash_seed_count,
                countmin_accum, num_inputs, top_values, top_freqs);
        garbleIO->flush();

        clear_fmpz_array(bucket0, num_sh);
        clear_fmpz_array(bucket1, num_sh);

        std::cout << "Top K = " << K << " values and freqs, decreasing\n";
        for (unsigned int i = 0; i < K; i++) {
          std::cout << "Value: " << top_values[i] << ", freq: " << top_freqs[i] << std::endl;
        }

        std::cout << "eval time: " << sec_from(start2) << std::endl;

        delete[] top_values;
        delete[] top_freqs;

        return RET_ANS;

    } else {
        size_t num_inputs;
        recv_size(serverfd, num_inputs);
        std::cout << "num_inputs: " << num_inputs << std::endl;

        std::string* const tag_list = new std::string[num_inputs];
        bool* const valid = new bool[num_inputs];
        bool* const shares_sh_x = new bool[num_inputs * share_size_sh];
        bool* const shares_sh_y = new bool[num_inputs * share_size_sh];
        bool* const shares_bucket = new bool[num_inputs * share_size_bucket];
        bool* const shares_layer = new bool[num_inputs * share_size_layer];
        bool* const shares_count = new bool[num_inputs * share_size_count];

        // Align by tag
        for (unsigned int i = 0; i < num_inputs; i++) {
            const std::string tag = get_tag(serverfd);
            tag_list[i] = tag;
            valid[i] = (share_map.find(tag) != share_map.end());
        }
        send_bool_batch(serverfd, valid, num_inputs);
        std::cout << "tag time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        for (unsigned int i = 0; i < num_inputs; i++) {
            if (!valid[i]) continue;
            bool* a; bool* b; bool* c; bool* d; bool* e;
            std::tie(a, b, c, d, e) = share_map[tag_list[i]];
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
        delete[] tag_list;

        auto start3 = clock_start();
        // Round 1: All B2A: (count, bucket, y)
        // count used (and freq checked), so first
        // bucket similarly just freq checked, so second
        // y final, with offset
        const size_t convert_size = num_inputs * share_convert_size;
        fmpz_t* shares_p; new_fmpz_array(&shares_p, convert_size);
        bool* shares_2 = new bool[convert_size];
        // Count first
        memcpy(shares_2, shares_count, num_inputs * share_size_count);
        delete[] shares_count;
        // Bucket second
        memcpy(&shares_2[num_inputs * share_size_count],
                shares_bucket, num_inputs * share_size_bucket);
        // Y third
        const size_t share_y_offset = num_inputs * (share_size_count + share_size_bucket);
        memcpy(&shares_2[share_y_offset], shares_sh_y, num_inputs * share_size_sh);
        delete[] shares_sh_y;
        sent_bytes += correlated_store->b2a_single(convert_size, shares_2, shares_p);
        std::cout << "B2A time: " << sec_from(start3) << std::endl; start3 = clock_start();

        // Round 1.5: Mask multiply (AND)
        // TODO: fold into above
        bool* mask = new bool[num_inputs * share_size_bucket * share_size_layer];
        sent_bytes += correlated_store->multiply_BoolShares_cross(
                num_inputs, share_size_bucket, share_size_layer,
                shares_bucket, shares_layer, mask);
        std::cout << "B2A time: " << sec_from(start3) << std::endl; start3 = clock_start();

        // Round 2: Frequency checks
        bool* parity = new bool[cfg.Q + cfg.countmin_cfg.d];
        fmpz_t* sums; new_fmpz_array(&sums, cfg.Q + cfg.countmin_cfg.d);
        size_t idx = 0;
        for (unsigned int i = 0; i < num_inputs; i++) {
            for (unsigned int j = 0; j < cfg.countmin_cfg.d; j++) {
                for (unsigned int k = 0; k < cfg.countmin_cfg.w; k++) {
                    parity[j + cfg.Q] ^= shares_2[idx];
                    fmpz_mod_add(sums[j], sums[j], shares_p[idx], mod_ctx);
                    idx++;
                }
            }
        }
        for (unsigned int i = 0; i < num_inputs; i++) {
            for (unsigned int j = 0; j < cfg.Q; j++) {
                for (unsigned int k = 0; k < cfg.B; k++) {
                    parity[j] ^= shares_2[idx];
                    fmpz_mod_add(sums[j + cfg.countmin_cfg.d], sums[j+ cfg.countmin_cfg.d], shares_p[idx], mod_ctx);
                    idx++;
                }
            }
        }
        delete[] shares_2;
        bool* parity_other = new bool[cfg.Q + cfg.countmin_cfg.d];
        fmpz_t* sums_other; new_fmpz_array(&sums_other, cfg.Q + cfg.countmin_cfg.d);


        sent_bytes += send_bool_batch(serverfd, parity, cfg.Q + cfg.countmin_cfg.d);
        sent_bytes += send_fmpz_batch(serverfd, sums, cfg.Q + cfg.countmin_cfg.d);
        recv_bool_batch(serverfd, parity_other, cfg.Q + cfg.countmin_cfg.d);
        recv_fmpz_batch(serverfd, sums_other, cfg.Q + cfg.countmin_cfg.d);

        bool all_valid = false;
        for (unsigned int j = 0; j < cfg.Q + cfg.countmin_cfg.d; j++) {
            all_valid |= (parity[j] ^ parity_other[j]) == (total_inputs % 2);
            fmpz_mod_add(sums[j], sums[j], sums_other[j], mod_ctx);
            all_valid |= fmpz_equal_ui(sums[j], total_inputs);
            if (!all_valid)
                continue;
        }
        delete[] parity;
        delete[] parity_other;
        clear_fmpz_array(sums, cfg.Q + cfg.countmin_cfg.d);
        clear_fmpz_array(sums_other, cfg.Q + cfg.countmin_cfg.d);
        if (all_valid) {
            memset(valid, true, num_inputs);
        } else {
            memset(valid, false, num_inputs);
            std::cout << "Batch not valid. Individual check currently not implemented" << std::endl;
        }
        std::cout << "freq validation time: " << sec_from(start3) << std::endl;

        // Round 2.5: Heavy convert with mask
        // 2 rounds of OT, which can't easily fold
        // Also accumulates into buckets
        start3 = clock_start();
        sent_bytes += correlated_store->heavy_convert_mask(
            num_inputs, cfg.Q, cfg.B * cfg.R, cfg.SH_depth,
            shares_sh_x, &shares_p[share_y_offset], mask, valid, bucket0, bucket1);
        // Also count-min accumulation
        fmpz_t* countmin_accum; new_fmpz_array(&countmin_accum, share_size_count);
        accumulate(num_inputs, share_size_count, shares_p, valid, countmin_accum);
        clear_fmpz_array(shares_p, convert_size);
        delete[] shares_sh_x;
        delete[] mask;
        std::cout << "heavy_convert time: " << sec_from(start3) << std::endl;
        std::cout << "convert+accum time: " << sec_from(start2) << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;
        std::cout << "compute bytes sent: " << sent_bytes << std::endl;

        // Evaluation: Validity check
        size_t num_valid = 0;
        for (unsigned int i = 0; i < num_inputs; i++)
            num_valid += valid[i];
        std::cout << "Final valid count: " << num_valid << " / " << total_inputs << std::endl;
        delete[] valid;
        if (num_valid < total_inputs * (1 - INVALID_THRESHOLD)) {
            std::cout << "Failing, This is less than the invalid threshold of " << INVALID_THRESHOLD << std::endl;
            return RET_INVALID;
        }

        // Evaluation: Answer extraction

        start2 = clock_start();

        uint64_t* top_values = new uint64_t[K];
        uint64_t* top_freqs = new uint64_t[K];

        full_heavy_extract(server_num, cfg, bucket0, bucket1, hash_seed_split, hash_seed_count,
                countmin_accum, num_inputs, top_values, top_freqs);
        garbleIO->flush();

        clear_fmpz_array(bucket0, num_sh);
        clear_fmpz_array(bucket1, num_sh);

        std::cout << "Top K = " << K << " values and freqs, decreasing\n";
        for (unsigned int i = 0; i < K; i++) {
          std::cout << "Value: " << top_values[i] << ", freq: " << top_freqs[i] << std::endl;
        }

        delete[] top_values;
        delete[] top_freqs;

        std::cout << "eval time: " << sec_from(start2) << std::endl;

        return RET_ANS;
    }
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
        correlated_store = new PrecomputeStore(serverfd, server_num, ot0, ot1, CACHE_SIZE, LAZY_PRECOMPUTE);
    // } else if (STORE_TYPE == ot_store) {
    //     correlated_store = new OTCorrelatedStore(serverfd, server_num, ot0, ot1);
    } else if (STORE_TYPE == validate_store) {
        correlated_store = new ValidateCorrelatedStore(serverfd, server_num, ot0, ot1, CACHE_SIZE, LAZY_PRECOMPUTE);
    } else {
        error_exit("Unknown/unsupported store type");
    }

    // Reuse IO from OT
    garbleIO = ot0->io;
    setup_semi_honest(garbleIO, server_num + 1);

    int sockfd, newsockfd;
    sockaddr_in addr;


    bind_and_listen(addr, sockfd, client_port, 1);

    while(1) {
        if (STORE_TYPE != ot_store) {
            ((PrecomputeStore*) correlated_store)->maybe_Update();
            ((PrecomputeStore*) correlated_store)->print_Sizes();
        }

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
