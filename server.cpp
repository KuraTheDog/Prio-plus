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

uint64_t randx_uses = 0;
fmpz_t randomX;
// Precomputes for the current random X, keyed by number of mults.
std::unordered_map<size_t, MultCheckPreComp*> eval_precomp_store;

OT_Wrapper* ot0;
OT_Wrapper* ot1;

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

// TODO: Consider marking all "invalid" paths as [[unlikely]], mainly for loops

size_t send_out(const int sockfd, const void* const buf, const size_t len) {
    size_t ret = send(sockfd, buf, len, 0);
    if (ret <= 0) error_exit("Failed to send");
    return ret;
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

// Gets a new randomX, syncs it between servers.
void sync_randomX(const int serverfd, const int server_num, fmpz_t randomX) {
    if (server_num == 0) {
        recv_fmpz(serverfd, randomX);
        // std::cout << "Got randomX: "; fmpz_print(randomX); std::cout << std::endl;
    } else {
        fmpz_randm(randomX, seed, Int_Modulus);
        // std::cout << "Sending randomX: "; fmpz_print(randomX); std::cout << std::endl;
        send_fmpz(serverfd, randomX);
    }
}

// TODO: can maybe batch this? I.e. get list of all PK at once.
std::string get_pk(const int serverfd) {
    char pk_buf[PK_LENGTH];
    recv_in(serverfd, &pk_buf[0], PK_LENGTH);
    std::string pk(pk_buf, pk_buf + PK_LENGTH);
    return pk;
}

MultCheckPreComp* getPrecomp(const size_t N) {
    MultCheckPreComp* pre;
    if (eval_precomp_store.find(N) == eval_precomp_store.end()) {
        pre = new MultCheckPreComp(N);
        pre->setEvalPoint(randomX);
        eval_precomp_store[N] = pre;
    } else {
        pre = eval_precomp_store[N];
    }
    return pre;
}

int recv_unvalidated(const int clientfd, const std::string pk, const size_t n) {
    if (STORE_TYPE != validate_store) {
        return 0;
    }

    // Overkill? Can be more efficient in the future. See client for more.
    const size_t N = NextPowerOfTwo(n-1);

    int num_bytes = 0;

    DaBit** const bits = new DaBit*[N];
    AltTriple** const trips = new AltTriple*[N];
    for (unsigned int i = 0; i < N; i++) {
        bits[i] = new DaBit();
        trips[i] = new AltTriple();
    }

    num_bytes += recv_DaBit_batch(clientfd, bits, N);
    num_bytes += recv_AltTriple_batch(clientfd, trips, N);

    ((ValidateCorrelatedStore*) correlated_store)->queueUnvalidated(bits, trips, pk);

    // std::cout << "got " << N << " unvalidated in " << num_bytes << " bytes" << std::endl;

    return num_bytes;
}

void process_unvalidated(const std::string pk, const size_t n) {
    if (STORE_TYPE != validate_store) {
        return;
    }
    const size_t N = NextPowerOfTwo(n-1);
    ((ValidateCorrelatedStore*) correlated_store)->processUnvalidated(pk, N);
}

// Currently shares_2 and shares_p are flat num_shares*num_values array.
// TODO: Consider reworking for matrix form
// TODO: return sent_bytes
fmpz_t* const share_convert(const size_t num_shares,  // # inputs
                            const size_t num_values,  // # values per input
                            const size_t* const num_bits,  // # bits per value
                            const uint64_t* const shares_2
                            ) {
    auto start = clock_start();
    [[maybe_unused]] int sent_bytes = 0;

    fmpz_t* shares_p; new_fmpz_array(&shares_p, num_shares * num_values);

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

    return shares_p;
}

// Batch of N (snips + num_input wire/share) validations
// Due to the nature of the final swap, both servers get the same valid array
const bool* const validate_snips(const size_t N,
                                 const size_t num_inputs,
                                 const int serverfd,
                                 const int server_num,
                                 Circuit* const * const circuit,
                                 const ClientPacket* const * const packet,
                                 const fmpz_t* const shares_p
                                 ) {
    auto start = clock_start();

    bool* const ans = new bool[N];

    const size_t NumRoots = NextPowerOfTwo(circuit[0]->NumMulGates());

    Checker** const checker = new Checker*[N];
    MultCheckPreComp* const pre = getPrecomp(NumRoots);
    randx_uses += N;
    for (unsigned int i = 0; i < N; i++)
        checker[i] = new Checker(circuit[i], server_num, packet[i], pre,
                                 &shares_p[i * num_inputs]);

    Cor** const cor = new Cor*[N];
    for (unsigned int i = 0; i < N; i++)
      cor[i] = checker[i]->CorFn();

    swap_Cor_batch(serverfd, cor, N);

    fmpz_t* valid_share; new_fmpz_array(&valid_share, N);
    for (unsigned int i = 0; i < N; i++) {
        checker[i]->OutShare(valid_share[i], cor[i]);
        delete cor[i];
        delete checker[i];
    }
    delete[] cor;
    delete[] checker;

    swap_fmpz_batch(serverfd, valid_share, N);

    for (unsigned int i = 0; i < N; i++) {
        ans[i] = fmpz_is_zero(valid_share[i]);
    }
    clear_fmpz_array(valid_share, N);

    std::cout << "snip circuit time: " << sec_from(start) << std::endl;
    return ans;
}

/*
shares_p: size inp X values, set across (all vals for inp0, then inp1, etc.)
valid: size inp, if false ignore corresponding shares
ans: size values, has modular sum of valid shares
returns number of valid.
*/
size_t accumulate(const size_t num_inputs,
                  const size_t num_values,
                  const fmpz_t* const shares_p,
                  const bool* const valid,
                  fmpz_t* const ans
                  ) {
    size_t num_valid = 0;

    for (unsigned int j = 0; j < num_values; j++)
        fmpz_zero(ans[j]);

    for (unsigned int i = 0; i < num_inputs; i++) {
        if (!valid[i])
            continue;
        for (unsigned int j = 0; j < num_values; j++) {
            fmpz_mod_add(ans[j], ans[j], shares_p[i * num_values + j], mod_ctx);
        }
        num_valid++;
    }

    return num_valid;
}

// Note: since bits, uses specific bitsum_ot, rather than normal store stuff.
returnType bit_sum(const initMsg msg, const int clientfd, const int serverfd,
                   const int server_num, uint64_t& ans) {
    std::unordered_map<std::string, bool> share_map;
    auto start = clock_start();

    BitShare share;
    const unsigned int total_inputs = msg.num_of_inputs;

    int num_bytes = 0;
    for (unsigned int i = 0; i < total_inputs; i++) {
        num_bytes += recv_in(clientfd, &share, sizeof(BitShare));
        const std::string pk(share.pk, share.pk + PK_LENGTH);
        // recv_unvalidated(clientfd, 1, pk);
        if (share_map.find(pk) != share_map.end())
            continue;
        share_map[pk] = share.val;
    }

    std::cout << "Received " << total_inputs << " total shares" << std::endl;
    std::cout << "bytes from client: " << num_bytes << std::endl;
    std::cout << "receive time: " << sec_from(start) << std::endl;
    start = clock_start();
    auto start2 = clock_start();

    // if (STORE_TYPE != ot_store)
    //     ((CorrelatedStore*) correlated_store)->checkDaBits(total_inputs);

    int server_bytes = 0;

    if (server_num == 1) {
        const size_t num_inputs = share_map.size();
        server_bytes += send_size(serverfd, num_inputs);
        bool* const shares = new bool[num_inputs];
        int i = 0;
        for (const auto& share : share_map) {
            server_bytes += send_out(serverfd, &share.first[0], PK_LENGTH);
            shares[i] = share.second;
            i++;
        }
        std::cout << "pk time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        const uint64_t b = bitsum_ot_receiver(ot0, shares, num_inputs);
        delete[] shares;

        send_uint64(serverfd, b);
        std::cout << "convert time: " << sec_from(start2) << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;
        std::cout << "sent server bytes: " << server_bytes << std::endl;
        return RET_NO_ANS;
    } else {
        size_t num_inputs, num_valid = 0;
        recv_size(serverfd, num_inputs);
        bool* const shares = new bool[num_inputs];
        bool* const valid = new bool[num_inputs];

        for (unsigned int i = 0; i < num_inputs; i++) {
            const std::string pk = get_pk(serverfd);

            bool is_valid = (share_map.find(pk) != share_map.end());
            valid[i] = is_valid;
            if (!is_valid)
                continue;
            num_valid++;
            shares[i] = share_map[pk];
        }
        std::cout << "pk time: " << sec_from(start2) << std::endl;
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

    int num_bytes = 0;
    for (unsigned int i = 0; i < total_inputs; i++) {
        num_bytes += recv_in(clientfd, &share, sizeof(IntShare));
        const std::string pk(share.pk, share.pk + PK_LENGTH);

        if (share_map.find(pk) != share_map.end()
            or share.val >= max_val)
            continue;
        share_map[pk] = share.val;

        // std::cout << "share[" << i << "] = " << share.val << std::endl;

        num_bytes += recv_unvalidated(clientfd, pk, msg.num_bits);
    }

    std::cout << "Received " << total_inputs << " total shares" << std::endl;
    std::cout << "bytes from client: " << num_bytes << std::endl;
    std::cout << "receive time: " << sec_from(start) << std::endl;

    if (STORE_TYPE != ot_store)
        ((CorrelatedStore*) correlated_store)->checkDaBits(total_inputs * msg.num_bits);

    start = clock_start();
    auto start2 = clock_start();

    int server_bytes = 0;

    if (server_num == 1) {
        const size_t num_inputs = share_map.size();
        server_bytes += send_size(serverfd, num_inputs);
        uint64_t* const shares = new uint64_t[num_inputs];
        int i = 0;
        for (const auto& share : share_map) {
            server_bytes += send_out(serverfd, &share.first[0], PK_LENGTH);
            shares[i] = share.second;
            i++;

            process_unvalidated(&share.first[0], msg.num_bits);
        }
        std::cout << "PK time: " << sec_from(start2) << std::endl;
        if (STORE_TYPE == validate_store) ((CorrelatedStore*) correlated_store)->checkDaBits(total_inputs * msg.num_bits);
        start2 = clock_start();
        fmpz_t* const shares_p = share_convert(num_inputs, 1, nbits, shares);
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
        std::cout << "sent server bytes: " << server_bytes << std::endl;
        return RET_NO_ANS;
    } else {
        size_t num_inputs;
        recv_size(serverfd, num_inputs);
        uint64_t* const shares = new uint64_t[num_inputs];
        bool* const valid = new bool[num_inputs];

        for (unsigned int i = 0; i < num_inputs; i++) {
            const std::string pk = get_pk(serverfd);

            bool is_valid = (share_map.find(pk) != share_map.end());
            valid[i] = is_valid;
            if (!is_valid)
                continue;
            shares[i] = share_map[pk];

            process_unvalidated(pk, msg.num_bits);
        }
        std::cout << "PK time: " << sec_from(start2) << std::endl;
        if (STORE_TYPE == validate_store) ((CorrelatedStore*) correlated_store)->checkDaBits(total_inputs * msg.num_bits);
        start2 = clock_start();
        fmpz_t* const shares_p = share_convert(num_inputs, 1, nbits, shares);
        std::cout << "convert time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        server_bytes += send_bool_batch(serverfd, valid, num_inputs);

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
        std::cout << "sent server bytes: " << server_bytes << std::endl;
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

    int num_bytes = 0;
    for (unsigned int i = 0; i < total_inputs; i++) {
        num_bytes += recv_in(clientfd, &share, sizeof(IntShare));
        const std::string pk(share.pk, share.pk + PK_LENGTH);

        if (share_map.find(pk) != share_map.end())
            continue;
        share_map[pk] = share.val;
    }

    std::cout << "Received " << total_inputs << " total shares" << std::endl;
    std::cout << "bytes from client: " << num_bytes << std::endl;
    std::cout << "receive time: " << sec_from(start) << std::endl;
    start = clock_start();
    auto start2 = clock_start();

    int server_bytes = 0;

    if (server_num == 1) {
        const size_t num_inputs = share_map.size();
        server_bytes += send_size(serverfd, num_inputs);
        uint64_t b = 0;
        std::string* const pk_list = new std::string[num_inputs];
        size_t idx = 0;
        for (const auto& share : share_map) {
            server_bytes += send_out(serverfd, &share.first[0], PK_LENGTH);
            pk_list[idx] = share.first;
            idx++;
        }
        std::cout << "PK time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        bool* const other_valid = new bool[num_inputs];
        recv_bool_batch(serverfd, other_valid, num_inputs);
        for (unsigned int i = 0; i < num_inputs; i++) {
            if (!other_valid[i])
                continue;
            b ^= share_map[pk_list[i]];
        }
        delete[] other_valid;

        send_uint64(serverfd, b);
        delete[] pk_list;
        std::cout << "convert time: " << sec_from(start2) << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;
        std::cout << "sent server bytes: " << server_bytes << std::endl;
        return RET_NO_ANS;
    } else {
        size_t num_inputs, num_valid = 0;
        recv_size(serverfd, num_inputs);
        uint64_t a = 0;
        bool* const valid = new bool[num_inputs];

        for (unsigned int i = 0; i < num_inputs; i++) {
            const std::string pk = get_pk(serverfd);
            valid[i] = (share_map.find(pk) != share_map.end());
            if (!valid[i])
                continue;
            num_valid++;
            a ^= share_map[pk];
        }

        std::cout << "PK + convert time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        server_bytes += send_bool_batch(serverfd, valid, num_inputs);

        delete[] valid;

        uint64_t b;
        recv_uint64(serverfd, b);
        const uint64_t aggr = a ^ b;

        std::cout << "Final valid count: " << num_valid << " / " << total_inputs << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;
        std::cout << "sent server bytes: " << server_bytes << std::endl;
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

    int num_bytes = 0;
    for (unsigned int i = 0; i < total_inputs; i++) {
        num_bytes += recv_in(clientfd, &share.pk[0], PK_LENGTH);
        const std::string pk(share.pk, share.pk + PK_LENGTH);

        num_bytes += recv_uint64_batch(clientfd, &shares[i*(B+1)], B+1);

        if (share_map.find(pk) != share_map.end())
            continue;
        share_map[pk] = &shares[i*(B+1)];
    }

    std::cout << "Received " << total_inputs << " total shares" << std::endl;
    std::cout << "bytes from client: " << num_bytes << std::endl;
    std::cout << "receive time: " << sec_from(start) << std::endl;
    start = clock_start();
    auto start2 = clock_start();

    int server_bytes = 0;

    if (server_num == 1) {
        const size_t num_inputs = share_map.size();
        server_bytes += send_size(serverfd, num_inputs);
        uint64_t b[B+1];
        memset(b, 0, sizeof(b));
        std::string* const pk_list = new std::string[num_inputs];
        size_t idx = 0;
        for (const auto& share : share_map) {
            server_bytes += send_out(serverfd, &share.first[0], PK_LENGTH);
            pk_list[idx] = share.first;
            idx++;
        }
        std::cout << "PK time: " << sec_from(start2) << std::endl;
        start2 = clock_start();
        bool* const other_valid = new bool[num_inputs];
        recv_bool_batch(serverfd, other_valid, num_inputs);
        for (unsigned int i = 0; i < num_inputs; i++) {
            if (!other_valid[i])
                continue;
            for (unsigned int j = 0; j <= B; j++)
                b[j] ^= share_map[pk_list[i]][j];
        }
        delete[] other_valid;
        delete[] shares;
        delete[] pk_list;
        send_uint64_batch(serverfd, b, B+1);

        std::cout << "convert time: " << sec_from(start2) << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;
        std::cout << "sent server bytes: " << server_bytes << std::endl;
        return RET_NO_ANS;
    } else {
        size_t num_inputs, num_valid = 0;
        recv_size(serverfd, num_inputs);
        uint64_t a[B+1];
        memset(a, 0, sizeof(a));
        bool* const valid = new bool[num_inputs];

        for (unsigned int i = 0; i < num_inputs; i++) {
            const std::string pk = get_pk(serverfd);
            valid[i] = (share_map.find(pk) != share_map.end());
            if (!valid[i])
                continue;
            num_valid++;
            for (unsigned int j = 0; j <= B; j++)
                a[j] ^= share_map[pk][j];
        }

        std::cout << "PK+convert time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        server_bytes += send_bool_batch(serverfd, valid, num_inputs);

        delete[] shares;
        delete[] valid;
        uint64_t b[B+1];
        recv_uint64_batch(serverfd, b, B+1);

        std::cout << "Final valid count: " << num_valid << " / " << total_inputs << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;
        std::cout << "sent server bytes: " << server_bytes << std::endl;
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

    int num_bytes = 0;
    for (unsigned int i = 0; i < total_inputs; i++) {
        num_bytes += recv_in(clientfd, &share, sizeof(VarShare));
        const std::string pk(share.pk, share.pk + PK_LENGTH);

        ClientPacket* const packet = new ClientPacket(NMul);
        int packet_bytes = recv_ClientPacket(clientfd, packet, NMul);
        num_bytes += packet_bytes;

        // std::cout << "share[" << i << "] = " << share.val << ", " << share.val_squared << std::endl;

        if ((share_map.find(pk) != share_map.end())
            or (share.val >= max_val)
            or (share.val_squared >= max_val * max_val)
            or (packet_bytes <= 0)
            ) {
            delete packet;
            continue;
        }
        share_map[pk] = {share.val, share.val_squared, packet};

        num_bytes += recv_unvalidated(clientfd, pk, msg.num_bits * num_dabits);
    }

    std::cout << "Received " << total_inputs << " total shares" << std::endl;
    std::cout << "bytes from client: " << num_bytes << std::endl;
    std::cout << "receive time: " << sec_from(start) << std::endl;

    if (STORE_TYPE != ot_store)
        ((CorrelatedStore*) correlated_store)->checkDaBits(total_inputs * msg.num_bits * num_dabits);

    start = clock_start();

    auto start2 = clock_start();

    int server_bytes = 0;

    if (server_num == 1) {
        const size_t num_inputs = share_map.size();
        server_bytes += send_size(serverfd, num_inputs);

        uint64_t* const shares = new uint64_t[2 * num_inputs];
        ClientPacket** const packet = new ClientPacket*[num_inputs];

        std::string* const pk_list = new std::string[num_inputs];
        Circuit** const circuit = new Circuit*[num_inputs];

        size_t idx = 0;
        for (const auto& share : share_map) {
            server_bytes += send_out(serverfd, &share.first[0], PK_LENGTH);
            pk_list[idx] = share.first;
            idx++;

            process_unvalidated(&share.first[0], msg.num_bits * num_dabits);
        }
        std::cout << "PK time: " << sec_from(start2) << std::endl;
        if (STORE_TYPE == validate_store) ((CorrelatedStore*) correlated_store)->checkDaBits(total_inputs * msg.num_bits * num_dabits);
        start2 = clock_start();

        for (unsigned int i = 0; i < num_inputs; i++) {
            uint64_t val = 0, val2 = 0;
            std::tie(val, val2, packet[i]) = share_map[pk_list[i]];
            circuit[i] = CheckVar();
            shares[2 * i] = val;
            shares[2 * i + 1] = val2;
        }
        fmpz_t* const shares_p = share_convert(num_inputs, 2,
                                               nbits, shares);
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
        delete[] pk_list;
        delete[] circuit;
        delete[] packet;
        delete[] shares;
        std::cout << "validate time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        // Convert
        fmpz_t* b; new_fmpz_array(&b, 2);
        accumulate(num_inputs, 2, shares_p, valid, b);

        std::cout << "accumulate time: " << sec_from(start2) << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;

        send_fmpz(serverfd, b[0]);
        send_fmpz(serverfd, b[1]);

        delete[] valid;
        clear_fmpz_array(b, 2);
        clear_fmpz_array(shares_p, num_inputs * 2);

        std::cout << "sent non-snip server bytes: " << server_bytes << std::endl;
        return RET_NO_ANS;
    } else {
        size_t num_inputs;
        recv_size(serverfd, num_inputs);

        uint64_t* const shares = new uint64_t[2 * num_inputs];
        ClientPacket** const packet = new ClientPacket*[num_inputs];

        bool* const valid = new bool[num_inputs];
        std::string* const pk_list = new std::string[num_inputs];
        Circuit** const circuit = new Circuit*[num_inputs];

        for (unsigned int i = 0; i < num_inputs; i++) {
            const std::string pk = get_pk(serverfd);
            pk_list[i] = pk;
            valid[i] = (share_map.find(pk) != share_map.end());

            process_unvalidated(pk, msg.num_bits * num_dabits);
        }
        std::cout << "PK time: " << sec_from(start2) << std::endl;
        if (STORE_TYPE == validate_store) ((CorrelatedStore*) correlated_store)->checkDaBits(total_inputs * msg.num_bits * num_dabits);
        start2 = clock_start();
        for (unsigned int i = 0; i < num_inputs; i++) {
            uint64_t val = 0, val2 = 0;
            if (valid[i]) {
                std::tie(val, val2, packet[i]) = share_map[pk_list[i]];
            } else {
                packet[i] = new ClientPacket(NMul);  // mock empty packet
            }
            circuit[i] = CheckVar();
            shares[2 * i] = val;
            shares[2 * i + 1] = val2;
        }
        fmpz_t* const shares_p = share_convert(num_inputs, 2,
                                               nbits, shares);
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
        server_bytes += send_bool_batch(serverfd, valid, num_inputs);
        delete[] snip_valid;
        delete[] pk_list;
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
        std::cout << "sent non-snip server bytes: " << server_bytes << std::endl;
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
    int num_bytes = 0;

    size_t degree;
    num_bytes += recv_size(clientfd, degree);

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

        num_bytes += recv_in(clientfd, &share.pk[0], PK_LENGTH);
        const std::string pk(share.pk, share.pk + PK_LENGTH);

        share.x_vals = new uint64_t[num_x];
        share.x2_vals = new uint64_t[num_quad];
        share.xy_vals = new uint64_t[num_x];

        num_bytes += recv_uint64_batch(clientfd, share.x_vals, num_x);
        num_bytes += recv_uint64(clientfd, share.y);
        num_bytes += recv_uint64_batch(clientfd, share.x2_vals, num_quad);
        num_bytes += recv_uint64_batch(clientfd, share.xy_vals, num_x);

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
        num_bytes += packet_bytes;

        if ((share_map.find(pk) != share_map.end())
            or (not sizes_valid)
            or (packet_bytes  <= 0)
            ) {
            delete[] share.x_vals;
            delete[] share.x2_vals;
            delete[] share.xy_vals;
            delete packet;
            continue;
        }

        share_map[pk] = {share.x_vals, share.y, share.x2_vals, share.xy_vals, packet};

        num_bytes += recv_unvalidated(clientfd, pk, msg.num_bits * num_dabits);
    }

    std::cout << "Received " << total_inputs << " total shares" << std::endl;
    std::cout << "bytes from client: " << num_bytes << std::endl;
    std::cout << "receive time: " << sec_from(start) << std::endl;

    if (STORE_TYPE != ot_store)
        ((CorrelatedStore*) correlated_store)->checkDaBits(total_inputs * msg.num_bits * num_dabits);

    start = clock_start();
    auto start2 = clock_start();

    int server_bytes = 0;

    if (server_num == 1) {
        const size_t num_inputs = share_map.size();
        server_bytes += send_size(serverfd, num_inputs);

        uint64_t* const shares = new uint64_t[num_inputs * num_fields];
        ClientPacket** const packet = new ClientPacket*[num_inputs];

        std::string* const pk_list = new std::string[num_inputs];
        Circuit** const circuit = new Circuit*[num_inputs];

        size_t idx = 0;
        for (const auto& share : share_map) {
            server_bytes += send_out(serverfd, &share.first[0], PK_LENGTH);
            pk_list[idx] = share.first;
            idx++;

            process_unvalidated(&share.first[0], msg.num_bits * num_dabits);
        }
        std::cout << "PK time: " << sec_from(start2) << std::endl;
        if (STORE_TYPE == validate_store) ((CorrelatedStore*) correlated_store)->checkDaBits(total_inputs * msg.num_bits * num_dabits);
        start2 = clock_start();

        for (unsigned int i = 0; i < num_inputs; i++) {
            uint64_t* x_vals;
            uint64_t y_val = 0;
            uint64_t* x2_vals;
            uint64_t* xy_vals;

            std::tie(x_vals, y_val, x2_vals, xy_vals, packet[i]) = share_map[pk_list[i]];
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
        fmpz_t* const shares_p = share_convert(num_inputs, num_fields,
                                               nbits, shares);
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
        delete[] pk_list;
        delete[] circuit;
        delete[] packet;
        delete[] shares;
        std::cout << "validate time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        // Convert
        fmpz_t* b; new_fmpz_array(&b, num_fields);
        accumulate(num_inputs, num_fields, shares_p, valid, b);

        std::cout << "accumulate time: " << sec_from(start2) << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;

        send_fmpz_batch(serverfd, b, num_fields);

        delete[] valid;
        clear_fmpz_array(b, num_fields);
        clear_fmpz_array(shares_p, num_inputs * num_fields);

        std::cout << "sent non-snip server bytes: " << server_bytes << std::endl;

        return RET_NO_ANS;
    } else {
        size_t num_inputs;
        recv_size(serverfd, num_inputs);

        uint64_t* const shares = new uint64_t[num_inputs * num_fields];
        ClientPacket** const packet = new ClientPacket*[num_inputs];

        bool* const valid = new bool[num_inputs];
        std::string* const pk_list = new std::string[num_inputs];
        Circuit** const circuit = new Circuit*[num_inputs];

        for (unsigned int i = 0; i < num_inputs; i++) {
            const std::string pk = get_pk(serverfd);
            pk_list[i] = pk;
            valid[i] = (share_map.find(pk) != share_map.end());

            process_unvalidated(pk, msg.num_bits * num_dabits);
        }
        std::cout << "PK time: " << sec_from(start2) << std::endl;
        if (STORE_TYPE == validate_store) ((CorrelatedStore*) correlated_store)->checkDaBits(total_inputs * msg.num_bits * num_dabits);
        start2 = clock_start();

        for (unsigned int i = 0; i < num_inputs; i++) {
            uint64_t* x_vals;
            uint64_t y_val = 0;
            uint64_t* x2_vals;
            uint64_t* xy_vals;
            if (valid[i]) {
                std::tie(x_vals, y_val, x2_vals, xy_vals, packet[i]) = share_map[pk_list[i]];
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
        fmpz_t* const shares_p = share_convert(num_inputs, num_fields,
                                               nbits, shares);
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
        server_bytes += send_bool_batch(serverfd, valid, num_inputs);
        delete[] snip_valid;
        delete[] pk_list;
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
        std::cout << "sent non-snip server bytes: " << server_bytes << std::endl;
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
    int num_bytes = 0;
    for (unsigned int i = 0; i < total_inputs; i++) {
        num_bytes += recv_in(clientfd, &share.pk[0], PK_LENGTH);
        const std::string pk(share.pk, share.pk+PK_LENGTH);
        share.arr = new bool[max_inp];
        num_bytes += recv_bool_batch(clientfd, share.arr, max_inp);

        if (share_map.find(pk) != share_map.end()) {
            delete[] share.arr;
            continue;
        }
        share_map[pk] = share.arr;

        // for (unsigned int j = 0; j < max_inp; j++) {
        //     std::cout << "share[" << i << ", " << j << "] = " << share.arr[j] << std::endl;
        // }
    }

    std::cout << "Received " << total_inputs << " total shares" << std::endl;
    std::cout << "bytes from client: " << num_bytes << std::endl;
    std::cout << "receive time: " << sec_from(start) << std::endl;

    if (STORE_TYPE != ot_store)
        ((CorrelatedStore*) correlated_store)->checkDaBits(total_inputs * max_inp);

    start = clock_start();
    auto start2 = clock_start();
    num_bytes = 0;

    if (server_num == 1) {
        const size_t num_inputs = share_map.size();
        num_bytes += send_size(serverfd, num_inputs);
        bool* const shares = new bool[num_inputs * max_inp];

        size_t idx = 0;
        for (const auto& share : share_map) {
            num_bytes += send_out(serverfd, &share.first[0], PK_LENGTH);
            memcpy(&shares[idx * max_inp], share.second, max_inp);
            delete[] share.second;
            idx++;
        }
        std::cout << "PK time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        fmpz_t* shares_p; new_fmpz_array(&shares_p, num_inputs * max_inp);
        num_bytes += correlated_store->b2a_single(num_inputs * max_inp, shares, shares_p);
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
        num_bytes += send_bool(serverfd, total_parity);
        num_bytes += send_fmpz(serverfd, sum);
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

        std::cout << "sent server bytes: " << num_bytes << std::endl;
        return RET_NO_ANS;
    } else {
        size_t num_inputs;
        recv_size(serverfd, num_inputs);
        bool* const shares = new bool[num_inputs * max_inp];
        bool* const valid = new bool[num_inputs];

        for (unsigned int i = 0; i < num_inputs; i++) {
            const std::string pk = get_pk(serverfd);
            valid[i] = (share_map.find(pk) != share_map.end());

            // realign shares_2 to pk order
            if (valid[i]) {
                memcpy(&shares[i * max_inp], share_map[pk], max_inp);
                delete[] share_map[pk];
            } else {
                memset(&shares[i * max_inp], 0, max_inp);
            }
        }
        std::cout << "PK time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        fmpz_t* shares_p; new_fmpz_array(&shares_p, num_inputs * max_inp);
        num_bytes += correlated_store->b2a_single(num_inputs * max_inp, shares, shares_p);
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
        num_bytes += send_bool(serverfd, total_parity);
        num_bytes += send_fmpz(serverfd, sum);
        recv_bool(serverfd, total_parity_other);
        recv_fmpz(serverfd, sum_other);
        bool all_valid = false;
        all_valid = (total_parity ^ total_parity_other) == (total_inputs % 2);
        if (all_valid) {
            fmpz_mod_add(sum, sum, sum_other, mod_ctx);
            all_valid = fmpz_equal_ui(sum, total_inputs);
        }
        num_bytes += send_bool(serverfd, all_valid);
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

            num_bytes += send_bool_batch(serverfd, valid, num_inputs);
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
        std::cout << "sent server bytes: " << num_bytes << std::endl;
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
    int num_bytes = 0;
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
        char pk_c[PK_LENGTH];
        num_bytes += recv_in(clientfd, &pk_c[0], PK_LENGTH);
        const std::string pk(pk_c, pk_c+PK_LENGTH);

        bool* const buff = new bool[2 * b];
        num_bytes += recv_bool_batch(clientfd, buff, 2 * b);

        if (share_map.find(pk) != share_map.end()) {
            delete[] buff;
            continue;
        }

        share_map[pk] = buff;
    }

    std::cout << "Received " << total_inputs << " total shares" << std::endl;
    std::cout << "bytes from client: " << num_bytes << std::endl;
    std::cout << "receive time: " << sec_from(start) << std::endl;

    correlated_store->checkDaBits(total_inputs * 2 * b);

    start = clock_start();
    auto start2 = clock_start();
    num_bytes = 0;
    if (server_num == 1) {
        const size_t num_inputs = share_map.size();
        const size_t n = num_inputs * b;
        num_bytes += send_size(serverfd, num_inputs);
        bool* const valid = new bool[num_inputs];
        bool* const x = new bool[n];
        bool* const y = new bool[n];

        size_t idx = 0;
        for (const auto& share : share_map) {
            num_bytes += send_out(serverfd, &share.first[0], PK_LENGTH);

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
        std::cout << "PK time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        num_bytes += correlated_store->heavy_convert(num_inputs, b, x, y, valid, bucket0, bucket1);
        delete[] x;
        delete[] y;
        std::cout << "convert+accum time: " << sec_from(start2) << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;
        std::cout << "compute bytes sent: " << num_bytes << std::endl;

        // Evaluate
        // True means |1| is larger than |0|
        start2 = clock_start();
        fmpz_t* larger; new_fmpz_array(&larger, b);
        num_bytes = correlated_store->abs_cmp(b, bucket0, bucket1, larger);
        num_bytes += send_fmpz_batch(serverfd, larger, b);
        std::cout << "evaluate time: " << sec_from(start2) << std::endl;
        std::cout << "evaluate bytes sent: " << num_bytes << std::endl;

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
            const std::string pk = get_pk(serverfd);
            valid[i] = (share_map.find(pk) != share_map.end());
            if (!valid[i]) {
                memset(&x[i * b], 0, b);
                memset(&y[i * b], 0, b);
                continue;
            }
            share = share_map[pk];
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
        num_bytes += send_bool_batch(serverfd, valid, num_inputs);
        std::cout << "PK time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        num_bytes += correlated_store->heavy_convert(num_inputs, b, x, y, valid, bucket0, bucket1);
        delete[] x;
        delete[] y;
        std::cout << "convert+accum time: " << sec_from(start2) << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;
        std::cout << "compute bytes sent: " << num_bytes << std::endl;

        // Evaluate

        start2 = clock_start();
        fmpz_t* larger_0; new_fmpz_array(&larger_0, b);
        fmpz_t* larger_1; new_fmpz_array(&larger_1, b);
        num_bytes = correlated_store->abs_cmp(b, bucket0, bucket1, larger_0);
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
        std::cout << "evaluate sent bytes: " << num_bytes << std::endl;

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

    int64_t num_bytes = 0;
    const size_t num_bits = msg.num_bits;
    const unsigned int total_inputs = msg.num_of_inputs;

    // Get other params
    size_t K; recv_size(clientfd, K);
    double delta; recv_double(clientfd, delta);
    double eps; recv_double(clientfd, eps);
    MultiHeavyConfig cfg(K, delta, num_bits, eps);
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
    const size_t num_sh = cfg.Q * cfg.B * cfg.SH_depth;
    // Which of the B SH instances it gets classified as
    HashStorePoly hash_classify(cfg.Q, num_bits, cfg.B, hash_seed_classify);
    // Split: each SH breakdown into the pairs, bucket 0 or 1. (original was by bits)
    // Base Q*B*depth, but can repeat across B. Invertible
    HashStoreBit hash_split(cfg.Q * cfg.SH_depth, num_bits, 2, hash_seed_split, cfg.SH_depth);
    // SingleHeavy +-1 values. 4-wize independent
    // Base Q*B*depth, but can repeat across B
    HashStorePoly hash_value(cfg.Q * cfg.SH_depth, num_bits, 2, hash_seed_value, 4);
    // SH storage: Q instances of B SH, each SH_depth large.
    fmpz_t* bucket0; new_fmpz_array(&bucket0, num_sh);
    fmpz_t* bucket1; new_fmpz_array(&bucket1, num_sh);

    for (unsigned int i = 0; i < total_inputs; i++) {
        char pk_c[PK_LENGTH];
        num_bytes += recv_in(clientfd, &pk_c[0], PK_LENGTH);
        const std::string pk(pk_c, pk_c+PK_LENGTH);

        bool* const share_sh_x = new bool[share_size_sh];
        bool* const share_sh_y = new bool[share_size_sh];
        bool* const share_mask = new bool[share_size_mask];
        bool* const share_count = new bool[share_size_count];
        num_bytes += recv_bool_batch(clientfd, share_sh_x, share_size_sh);
        num_bytes += recv_bool_batch(clientfd, share_sh_y, share_size_sh);
        num_bytes += recv_bool_batch(clientfd, share_mask, share_size_mask);
        num_bytes += recv_bool_batch(clientfd, share_count, share_size_count);

        if (share_map.find(pk) != share_map.end()) {
            delete[] share_sh_x;
            delete[] share_sh_y;
            delete[] share_mask;
            delete[] share_count;
            continue;
        }

        share_map[pk] = {share_sh_x, share_sh_y, share_mask, share_count};
    }

    std::cout << "Received " << total_inputs << " total shares" << std::endl;
    std::cout << "bytes from client: " << num_bytes << std::endl;
    std::cout << "receive time: " << sec_from(start) << std::endl;

    // For each pair of buckets, do 1 B2A.
    // All of count-min. Also all of mask for validation
    std::cout << "Inital dabit check" << std::endl;
    correlated_store->checkDaBits(total_inputs *
        (share_size_sh + share_size_count + share_size_mask));

    /* Stages:
    Validate each b (share_mask) is frequency vector. With B2A to check sum to 1
    "Validate H_q vector correct". share_sh. Not needed for R=1, is just entry
    (joint share validity)
    Update Countmin with share_count (and freq check in parallel)
    For each (q, b):
      final mask = share_sh & R_mask. Skip
      SH update, except OT multiply [z] by mask.

    mask: B2A just for validation. Kept as bool for Ot multiply
    sh: half B2A for use (z), half kept selection
    count: B2A for validation and accumulation
    */

    start = clock_start();
    auto start2 = clock_start();
    num_bytes = 0;
    if (server_num == 1) {
        const size_t num_inputs = share_map.size();
        num_bytes += send_size(serverfd, num_inputs);
        std::cout << "num_inputs: " << num_inputs << std::endl;

        std::string* const pk_list = new std::string[num_inputs];
        bool* const valid = new bool[num_inputs];
        bool* const shares_sh_x = new bool[num_inputs * share_size_sh];
        bool* const shares_sh_y = new bool[num_inputs * share_size_sh];
        bool* const shares_mask = new bool[num_inputs * share_size_mask];
        bool* const shares_count = new bool[num_inputs * share_size_count];

        size_t idx = 0;
        for (const auto& share : share_map) {
            num_bytes += send_out(serverfd, &share.first[0], PK_LENGTH);
            pk_list[idx] = share.first;
            idx++;
        }
        recv_bool_batch(serverfd, valid, num_inputs);
        std::cout << "PK time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        for (unsigned int i = 0; i < num_inputs; i++) {
            if (!valid[i]) continue;
            bool* a; bool* b; bool* c; bool* d;
            std::tie(a, b, c, d) = share_map[pk_list[i]];
            memcpy(&shares_sh_x[i * share_size_sh], a, share_size_sh);
            memcpy(&shares_sh_y[i * share_size_sh], b, share_size_sh);
            memcpy(&shares_mask[i * share_size_mask], c, share_size_mask);
            memcpy(&shares_count[i * share_size_count], d, share_size_count);
            delete a;
            delete b;
            delete c;
            delete d;
        }
        delete[] pk_list;

        auto start3 = clock_start();

        // Single round all B2A: count, mask, y
        const size_t convert_size = num_inputs * (
            share_size_count + share_size_mask + share_size_sh);
        fmpz_t* shares_p; new_fmpz_array(&shares_p, convert_size);
        bool* shares_2 = new bool[convert_size];
        memcpy(shares_2, shares_count, num_inputs * share_size_count);
        memcpy(&shares_2[num_inputs * (share_size_count)],
                shares_mask, num_inputs * share_size_mask);
        const size_t share_y_offset = num_inputs * (share_size_count + share_size_mask);
        memcpy(&shares_2[share_y_offset], shares_sh_y, num_inputs * share_size_sh);
        delete[] shares_sh_y;
        delete[] shares_count;
        num_bytes += correlated_store->b2a_single(convert_size, shares_2, shares_p);
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
        num_bytes += send_bool_batch(serverfd, parity, cfg.Q + cfg.countmin_cfg.d);
        num_bytes += send_fmpz_batch(serverfd, sums, cfg.Q + cfg.countmin_cfg.d);
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
        num_bytes += correlated_store->heavy_convert_mask(
            num_inputs, cfg.Q, cfg.B, cfg.SH_depth,
            shares_sh_x, &shares_p[share_y_offset], shares_mask,
            valid, bucket0, bucket1);
        clear_fmpz_array(shares_p, num_inputs * (share_size_count + share_size_mask));
        delete[] valid;
        delete[] shares_sh_x;
        delete[] shares_mask;
        std::cout << "heavy_convert time: " << sec_from(start3) << std::endl;
        std::cout << "convert+accum time: " << sec_from(start2) << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;
        std::cout << "compute bytes sent: " << num_bytes << std::endl;

        // straightforward eval. TODO: abs_cmp.
        start2 = clock_start();
        num_bytes += send_fmpz_batch(serverfd, bucket0, num_sh);
        num_bytes += send_fmpz_batch(serverfd, bucket1, num_sh);
        num_bytes += send_fmpz_batch(serverfd, countmin_accum, share_size_count);

        clear_fmpz_array(bucket0, num_sh);
        clear_fmpz_array(bucket1, num_sh);
        clear_fmpz_array(countmin_accum, share_size_count);

        return RET_NO_ANS;
    } else {
        size_t num_inputs;
        recv_size(serverfd, num_inputs);
        std::cout << "num_inputs: " << num_inputs << std::endl;

        std::string* const pk_list = new std::string[num_inputs];
        bool* const valid = new bool[num_inputs];
        bool* const shares_sh_x = new bool[num_inputs * share_size_sh];
        bool* const shares_sh_y = new bool[num_inputs * share_size_sh];
        bool* const shares_mask = new bool[num_inputs * share_size_mask];
        bool* const shares_count = new bool[num_inputs * share_size_count];

        for (unsigned int i = 0; i < num_inputs; i++) {
            const std::string pk = get_pk(serverfd);
            pk_list[i] = pk;
            valid[i] = (share_map.find(pk) != share_map.end());
        }
        send_bool_batch(serverfd, valid, num_inputs);
        std::cout << "PK time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        for (unsigned int i = 0; i < num_inputs; i++) {
            if (!valid[i]) continue;
            bool* a; bool* b; bool* c; bool* d;
            std::tie(a, b, c, d) = share_map[pk_list[i]];
            memcpy(&shares_sh_x[i * share_size_sh], a, share_size_sh);
            memcpy(&shares_sh_y[i * share_size_sh], b, share_size_sh);
            memcpy(&shares_mask[i * share_size_mask], c, share_size_mask);
            memcpy(&shares_count[i * share_size_count], d, share_size_count);
            delete a;
            delete b;
            delete c;
            delete d;
        }
        delete[] pk_list;

        auto start3 = clock_start();

        // Single round all B2A: count, mask, y
        const size_t convert_size = num_inputs * (
            share_size_count + share_size_mask + share_size_sh);
        fmpz_t* shares_p; new_fmpz_array(&shares_p, convert_size);
        bool* shares_2 = new bool[convert_size];
        memcpy(shares_2, shares_count, num_inputs * share_size_count);
        memcpy(&shares_2[num_inputs * (share_size_count)],
                shares_mask, num_inputs * share_size_mask);
        const size_t share_y_offset = num_inputs * (share_size_count + share_size_mask);
        memcpy(&shares_2[share_y_offset], shares_sh_y, num_inputs * share_size_sh);
        delete[] shares_sh_y;
        delete[] shares_count;
        num_bytes += correlated_store->b2a_single(convert_size, shares_2, shares_p);
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
        num_bytes += send_bool_batch(serverfd, parity, cfg.Q + cfg.countmin_cfg.d);
        num_bytes += send_fmpz_batch(serverfd, sums, cfg.Q + cfg.countmin_cfg.d);
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
        num_bytes += correlated_store->heavy_convert_mask(
            num_inputs, cfg.Q, cfg.B, cfg.SH_depth,
            shares_sh_x, &shares_p[share_y_offset], shares_mask,
            valid, bucket0, bucket1);
        clear_fmpz_array(shares_p, num_inputs * (share_size_count + share_size_mask));
        delete[] shares_sh_x;
        delete[] shares_mask;
        std::cout << "heavy_convert time: " << sec_from(start3) << std::endl;
        std::cout << "convert+accum time: " << sec_from(start2) << std::endl;
        std::cout << "total compute time: " << sec_from(start) << std::endl;
        std::cout << "compute bytes sent: " << num_bytes << std::endl;

        // Lazy eval. TODO: abs_cmp
        start2 = clock_start();

        fmpz_t* bucket0_other; new_fmpz_array(&bucket0_other, num_sh);
        fmpz_t* bucket1_other; new_fmpz_array(&bucket1_other, num_sh);
        fmpz_t* count_other; new_fmpz_array(&count_other, share_size_count);
        recv_fmpz_batch(serverfd, bucket0_other, num_sh);
        recv_fmpz_batch(serverfd, bucket1_other, num_sh);
        recv_fmpz_batch(serverfd, count_other, share_size_count);

        size_t num_valid = 0;
        for (unsigned int i = 0; i < num_inputs; i++)
            num_valid += valid[i];
        std::cout << "Final valid count: " << num_valid << " / " << total_inputs << std::endl;
        delete[] valid;
        if (num_valid < total_inputs * (1 - INVALID_THRESHOLD)) {
            std::cout << "Failing, This is less than the invalid threshold of " << INVALID_THRESHOLD << std::endl;
            return RET_INVALID;
        }

        fmpz_t* values; new_fmpz_array(&values, cfg.SH_depth);
        fmpz_t b0; fmpz_init(b0);
        fmpz_t b1; fmpz_init(b1);
        std::set<int64_t> candidates;
        for (unsigned int q = 0; q < cfg.Q; q++) {
            for (unsigned int b = 0; b < cfg.B; b++) {
                const size_t group_idx = q * cfg.B + b;
                for (unsigned int d = 0; d < cfg.SH_depth; d++) {
                    const size_t idx = group_idx * cfg.SH_depth + d;
                    fmpz_mod_add(b0, bucket0[idx], bucket0_other[idx], mod_ctx);
                    to_fsigned(b0, Int_Modulus);
                    fmpz_mod_add(b1, bucket1[idx], bucket1_other[idx], mod_ctx);
                    to_fsigned(b1, Int_Modulus);

                    // If b0 + b1 < some threshold, then ignore?
                    // Maybe e.g. b0 + b1 vs total / B (/ 2) ? (half of average?)
                    // Or some absolute thing
                    // Will be eliminated by large quantity, count-min, and validation anyways

                    fmpz_set_ui(values[d], fmpz_cmpabs(b0, b1) < 0 ? 1 : 0);
                    // std::cout << "Bucket[" << idx << "] (" << q << ", " << b << ", " << d << ") = ";
                    // std::cout << get_fsigned(b0, Int_Modulus) << ", " << get_fsigned(b1, Int_Modulus);
                    // std::cout << "  : cmp abs = " << fmpz_cmpabs(b0, b1);
                    // std::cout << ", value = " << fmpz_get_si(values[d]) << "\n";
                }
                uint64_t ans;
                int bad_hashes = hash_split.solve(q, values, ans);
                // std::cout << "Candidate[" << q << ", " << b << "] = ";
                // std::cout << ans << ", with " << bad_hashes << " invalid\n";
                if (bad_hashes == 0) {
                    candidates.insert(ans);
                }
            }
        }
        clear_fmpz_array(values, cfg.SH_depth);
        clear_fmpz_array(bucket0, num_sh);
        clear_fmpz_array(bucket1, num_sh);
        clear_fmpz_array(bucket0_other, num_sh);
        clear_fmpz_array(bucket1_other, num_sh);

        // Query count-min
        CountMin count_min(cfg.countmin_cfg);
        count_min.setStore(num_bits, hash_seed_count);
        count_min.init();

        for (unsigned int i = 0; i < share_size_count; i++) {
            fmpz_mod_add(count_min.counts[i], countmin_accum[i], count_other[i], mod_ctx);
        }
        clear_fmpz_array(countmin_accum, share_size_count);
        clear_fmpz_array(count_other, share_size_count);

        // std::cout << "combined countmin" << std::endl;
        // count_min.print();

        // Priority queue?
        // Size cap at total candidates, which is B * Q.
        // So don't need to worry about discarding too-small items as more come in
        std::priority_queue<std::pair<uint64_t, uint64_t>> frequencies;
        for (auto it = candidates.begin(); it!=candidates.end(); ++it) {
            uint64_t candidate = *it;
            uint64_t freq = count_min.query(candidate);
            frequencies.push(std::make_pair(freq, candidate));
        }

        std::cout << "Heavy items:" << std::endl;
        // Just print top 2K, for wiggle room etc.
        for (unsigned int i = 0; i < 2 * K; i++) {
        // while (!frequencies.empty()) {
            if (frequencies.empty()) break;
            auto next = frequencies.top();
            if (next.first == 0) break;
            std::cout << "\t value: " << next.second << "\t freq: " << next.first << std::endl;
            frequencies.pop();
        }

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

    fmpz_init(randomX);
    sync_randomX(serverfd, server_num, randomX);

    syncSnipSeeds(serverfd, server_num);

    ot0 = new OT_Wrapper(server_num == 0 ? nullptr : SERVER0_IP, 60051);
    ot1 = new OT_Wrapper(server_num == 1 ? nullptr : SERVER1_IP, 60052);

    // TODO: OT disabled for now.
    if (STORE_TYPE == precompute_store) {
        correlated_store = new PrecomputeStore(serverfd, server_num, ot0, ot1, CACHE_SIZE, LAZY_PRECOMPUTE);
    // } else if (STORE_TYPE == ot_store) {
    //     correlated_store = new OTCorrelatedStore(serverfd, server_num, ot0, ot1);
    } else if (STORE_TYPE == validate_store) {
        correlated_store = new ValidateCorrelatedStore(serverfd, server_num, ot0, ot1, CACHE_SIZE, LAZY_PRECOMPUTE);
    } else {
        error_exit("Unknown/unsupported store type");
    }

    int sockfd, newsockfd;
    sockaddr_in addr;

    bind_and_listen(addr, sockfd, client_port, 1);

    while(1) {
        // Refresh randomX if used too much
        if (randx_uses > EVAL_REUSE_THRESHOLD) {
            randx_uses = 0;
            sync_randomX(serverfd, server_num, randomX);
            // Update precomps
            for (const auto& pair : eval_precomp_store)
                pair.second -> setEvalPoint(randomX);
        }

        if (STORE_TYPE != ot_store) {
            ((PrecomputeStore*) correlated_store)->maybeUpdate();
            ((PrecomputeStore*) correlated_store)->printSizes();
        }

        socklen_t addrlen = sizeof(addr);

        std::cout << "waiting for connection..." << std::endl;

        newsockfd = accept(sockfd, (struct sockaddr*)&addr, &addrlen);
        if (newsockfd < 0) error_exit("Connection creation failure");

        // Get an initMsg
        initMsg msg;
        recv_in(newsockfd, &msg, sizeof(initMsg));

        if (msg.type == BIT_SUM) {
            std::cout << "BIT_SUM" << std::endl;
            auto start = clock_start();

            uint64_t ans;
            returnType ret = bit_sum(msg, newsockfd, serverfd, server_num, ans);
            if (ret == RET_ANS)
                std::cout << "Ans: " << ans << std::endl;

            std::cout << "Total time  : " << sec_from(start) << std::endl;
        } else if (msg.type == INT_SUM) {
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
        } else if (msg.type == NONE_OP) {
            std::cout << "Empty client message" << std::endl;
        } else {
            std::cout << "Unrecognized message type: " << msg.type << std::endl;
        }
        close(newsockfd);
    }

    delete correlated_store;
    for (const auto& precomp : eval_precomp_store)
        delete precomp.second;

    delete ot0;
    delete ot1;
    fmpz_clear(randomX);

    RootManager(1).clearCache();
    clear_constants();

    return 0;
}
