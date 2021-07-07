#include "server.h"

#include <netinet/in.h>
#include <sys/socket.h>
#include <sys/wait.h>
#include <unistd.h>

#include <cstdlib>
#include <iostream>
#include <unordered_map>
#include <string>

#include "correlated.h"
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

// Can keep the same random X for a while
#define RANDOMX_THRESHOLD 1e6
uint64_t randx_uses = 0;
fmpz_t randomX;
// Precomputes for the current random X
std::unordered_map<size_t, CheckerPreComp*> precomp_store;

OT_Wrapper* ot0;
OT_Wrapper* ot1;

// Precompute cache of edabits and beaver triples
CorrelatedStore* correlated_store;
// #define CACHE_SIZE 8192
// #define CACHE_SIZE 65536
// #define CACHE_SIZE 262144
#define CACHE_SIZE 2097152
// If set, does fast but insecure offline precompute.
#define LAZY_PRECOMPUTE true
// Generate excess in case it runs out
#define OVER_PRECOMPUTE false
// Whether to use OT or Dabits
#define USE_OT_B2A true

uint64_t int_sum_max;
uint32_t num_bits;

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

// TODO: can maybe batch this? Seems to help in other places. Not a big timesink though, so should be fine for now
std::string get_pk(const int serverfd) {
    char pk_buf[PK_LENGTH];
    recv_in(serverfd, &pk_buf[0], PK_LENGTH);
    std::string pk(pk_buf, pk_buf + PK_LENGTH);
    return pk;
}

CheckerPreComp* getPrecomp(const size_t N) {
    CheckerPreComp* pre;
    if (precomp_store.find(N) == precomp_store.end()) {
        pre = new CheckerPreComp(N);
        pre->setCheckerPrecomp(randomX);
        precomp_store[N] = pre;
    } else {
        pre = precomp_store[N];
    }
    return pre;
}

// Currently shares_2 and shares_p are flat num_shares*num_values array.
// TODO: Consider reworking for matrix form
fmpz_t* share_convert(const size_t num_shares,
                      const size_t num_values,
                      const size_t* const num_bits,
                      const uint64_t* const shares_2
                      ) {
    auto start = clock_start();

    fmpz_t* shares_p;

    // convert
    fmpz_t* f_shares2; new_fmpz_array(&f_shares2, num_shares * num_values);
    for (unsigned int i = 0; i < num_shares; i++) {
        for (unsigned int j = 0; j < num_values; j++)
            fmpz_set_ui(f_shares2[i * num_values + j],
                        shares_2[i * num_values + j]);
    }

    if (USE_OT_B2A) {
        const size_t mod = fmpz_get_ui(Int_Modulus);
        // TODO: maybe make it take not flattened?
        shares_p = correlated_store->b2a_ot(
            num_shares, num_values, num_bits, f_shares2, mod);

    } else {  // edabit conversion
        size_t* const bits_arr = new size_t[num_shares * num_values];
        for (unsigned int i = 0; i < num_shares; i++)
            memcpy(&bits_arr[i * num_values], num_bits, num_values * sizeof(size_t));

        shares_p = correlated_store->b2a_daBit_multi(
            num_shares * num_values, bits_arr, f_shares2);

        delete[] bits_arr;
    }

    clear_fmpz_array(f_shares2, num_shares * num_values);

    std::cout << "Share convert time: " << sec_from(start) << std::endl;

    return shares_p;
}

// Batch of N (snips + num_input wire/share) validations
// Due to the nature of the final swap, both servers get the same valid array
bool* validate_snips(const size_t N,
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
    pid_t pid = 0;
    int status = 0;

    init_roots(NumRoots);

    Checker** const checker = new Checker*[N];
    CheckerPreComp* const pre = getPrecomp(NumRoots);
    randx_uses += N;
    for (unsigned int i = 0; i < N; i++)
        checker[i] = new Checker(circuit[i], server_num, packet[i], pre,
                                 &shares_p[i * num_inputs]);

    CorShare** const cor_share = new CorShare*[N];
    for (unsigned int i = 0; i < N; i++)
      cor_share[i] = checker[i]->CorShareFn();

    if (correlated_store->do_fork) pid = fork();
    if (pid == 0) {
        send_CorShare_batch(serverfd, cor_share, N);
        if (correlated_store->do_fork) exit(EXIT_SUCCESS);
    }
    CorShare** const cor_share_other = new CorShare*[N];
    for (unsigned int i = 0; i < N; i++)
        cor_share_other[i] = new CorShare();
    recv_CorShare_batch(serverfd, cor_share_other, N);

    fmpz_t* valid_share; new_fmpz_array(&valid_share, N);
    for (unsigned int i = 0; i < N; i++) {
        const Cor* const cor = new Cor(cor_share[i], cor_share_other[i]);
        checker[i]->OutShare(valid_share[i], cor);
        delete cor;
    }
    if (correlated_store->do_fork) waitpid(pid, &status, 0);

    // TODO: Can be simplified: one sends share, other sends if valid
    if (correlated_store->do_fork) pid = fork();
    if (pid == 0) {
        send_fmpz_batch(serverfd, valid_share, N);
        if (correlated_store->do_fork) exit(EXIT_SUCCESS);
    }
    fmpz_t* valid_share_other; new_fmpz_array(&valid_share_other, N);
    recv_fmpz_batch(serverfd, valid_share_other, N);

    for (unsigned int i = 0; i < N; i++) {
        ans[i] = AddToZero(valid_share[i], valid_share_other[i]);
    }
    clear_fmpz_array(valid_share, N);
    clear_fmpz_array(valid_share_other, N);

    for (unsigned int i = 0; i < N; i++) {
        delete cor_share[i];
        delete cor_share_other[i];
        delete checker[i];
    }
    delete[] cor_share;
    delete[] cor_share_other;
    delete[] checker;

    if (correlated_store->do_fork) waitpid(pid, &status, 0);

    std::cout << "snip circuit time: " << sec_from(start) << std::endl;
    return ans;
}

fmpz_t* accumulate(const size_t N,
                   const size_t num_inputs,
                   const fmpz_t* const shares_p,
                   const bool* const valid
                   ) {
    fmpz_t* ans; new_fmpz_array(&ans, num_inputs);

    for (unsigned int j = 0; j < num_inputs; j++)
        fmpz_zero(ans[j]);

    for (unsigned int i = 0; i < N; i++) {
        if (!valid[i])
            continue;
        for (unsigned int j = 0; j < num_inputs; j++) {
            fmpz_add(ans[j], ans[j], shares_p[i * num_inputs + j]);
            fmpz_mod(ans[j], ans[j], Int_Modulus);  // How frequently?
        }
    }

    return ans;
}

returnType bit_sum(const initMsg msg, const int clientfd, const int serverfd, const int server_num, uint64_t& ans) {
    std::unordered_map<std::string, bool> share_map;
    auto start = clock_start();

    BitShare share;
    const unsigned int total_inputs = msg.num_of_inputs;

    int num_bytes = 0;
    for (unsigned int i = 0; i < total_inputs; i++) {
        num_bytes += recv_in(clientfd, &share, sizeof(BitShare));
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
        const unsigned int num_inputs = share_map.size();
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
        std::cout << "compute time: " << sec_from(start) << std::endl;
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
        std::cout << "compute time: " << sec_from(start) << std::endl;
        if (num_valid < total_inputs * (1 - INVALID_THRESHOLD)) {
            std::cout << "Failing, This is less than the invalid threshold of " << INVALID_THRESHOLD << std::endl;
            return RET_INVALID;
        }

        ans = a + b;
        return RET_ANS;
    }
}

returnType int_sum(const initMsg msg, const int clientfd, const int serverfd, const int server_num, uint64_t& ans) {
    std::unordered_map<std::string, uint64_t> share_map;
    auto start = clock_start();

    IntShare share;
    const uint64_t max_val = 1ULL << num_bits;
    const unsigned int total_inputs = msg.num_of_inputs;
    const size_t nbits[1] = {num_bits};

    int num_bytes = 0;
    for (unsigned int i = 0; i < total_inputs; i++) {
        num_bytes += recv_in(clientfd, &share, sizeof(IntShare));
        const std::string pk(share.pk, share.pk + PK_LENGTH);

        if (share_map.find(pk) != share_map.end()
            or share.val >= max_val)
            continue;
        share_map[pk] = share.val;

        // std::cout << "share[" << i << "] = " << share.val << std::endl;
    }

    std::cout << "Received " << total_inputs << " total shares" << std::endl;
    std::cout << "bytes from client: " << num_bytes << std::endl;
    std::cout << "receive time: " << sec_from(start) << std::endl;
    start = clock_start();
    auto start2 = clock_start();

    int server_bytes = 0;

    if (server_num == 1) {
        const unsigned int num_inputs = share_map.size();
        server_bytes += send_size(serverfd, num_inputs);
        uint64_t** const shares = new uint64_t*[num_inputs];
        int i = 0;
        for (const auto& share : share_map) {
            server_bytes += send_out(serverfd, &share.first[0], PK_LENGTH);
            shares[i] = new uint64_t[1];
            shares[i][0] = share.second;
            i++;
        }
        std::cout << "PK time: " << sec_from(start2) << std::endl;
        start2 = clock_start();
        const uint64_t* const * const b_all = intsum_ot_receiver(ot0, shares, nbits, num_inputs, 1);
        uint64_t b = 0;
        for (unsigned int i = 0; i < num_inputs; i++) {
            b += b_all[i][0];
            delete[] shares[i];
            delete[] b_all[i];
        }
        delete[] shares;
        delete[] b_all;

        send_uint64(serverfd, b);
        std::cout << "convert time: " << sec_from(start2) << std::endl;
        std::cout << "compute time: " << sec_from(start) << std::endl;
        std::cout << "sent server bytes: " << server_bytes << std::endl;
        return RET_NO_ANS;
    } else {
        size_t num_inputs, num_valid = 0;
        recv_size(serverfd, num_inputs);
        uint64_t** const shares = new uint64_t*[num_inputs];
        bool* const valid = new bool[num_inputs];

        for (unsigned int i = 0; i < num_inputs; i++) {
            const std::string pk = get_pk(serverfd);

            bool is_valid = (share_map.find(pk) != share_map.end());
            valid[i] = is_valid;
            shares[i] = new uint64_t[1];
            if (!is_valid)
                continue;
            num_valid++;
            shares[i][0] = share_map[pk];
        }
        std::cout << "PK time: " << sec_from(start2) << std::endl;
        start2 = clock_start();
        const uint64_t* const * const a_all = intsum_ot_sender(ot0, shares, valid, nbits, num_inputs, 1);
        delete[] valid;

        uint64_t a = 0;
        for (unsigned int i = 0; i < num_inputs; i++) {
            a += a_all[i][0];
            delete[] shares[i];
            delete[] a_all[i];
        }
        delete[] shares;
        delete[] a_all;

        uint64_t b;
        recv_uint64(serverfd, b);
        std::cout << "Final valid count: " << num_valid << " / " << total_inputs << std::endl;
        std::cout << "convert time: " << sec_from(start2) << std::endl;
        std::cout << "compute time: " << sec_from(start) << std::endl;
        if (num_valid < total_inputs * (1 - INVALID_THRESHOLD)) {
            std::cout << "Failing, This is less than the invalid threshold of " << INVALID_THRESHOLD << std::endl;
            return RET_INVALID;
        }

        ans = a + b;
        return RET_ANS;
    }
}

// For AND and OR
returnType xor_op(const initMsg msg, const int clientfd, const int serverfd, const int server_num, bool& ans) {
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
        const unsigned int num_inputs = share_map.size();
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
        std::cout << "compute time: " << sec_from(start) << std::endl;
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
        std::cout << "compute time: " << sec_from(start) << std::endl;
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
returnType max_op(const initMsg msg, const int clientfd, const int serverfd, const int server_num, uint64_t& ans) {
    std::unordered_map<std::string, uint64_t*> share_map;
    auto start = clock_start();

    MaxShare share;
    const unsigned int total_inputs = msg.num_of_inputs;
    const unsigned int B = msg.max_inp;
    const size_t share_sz = (B+1) * sizeof(uint64_t);
    // Need this to have all share arrays stay in memory, for server1 later.
    uint64_t* const shares = new uint64_t[total_inputs * (B + 1)];

    int num_bytes = 0;
    for (unsigned int i = 0; i < total_inputs; i++) {
        num_bytes += recv_in(clientfd, &share.pk[0], PK_LENGTH);
        const std::string pk(share.pk, share.pk + PK_LENGTH);

        num_bytes += recv_in(clientfd, &shares[i*(B+1)], share_sz);

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
        const unsigned int num_inputs = share_map.size();
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
        send_out(serverfd, &b[0], share_sz);
        delete[] shares;
        delete[] pk_list;

        std::cout << "convert time: " << sec_from(start2) << std::endl;
        std::cout << "compute time: " << sec_from(start) << std::endl;
        std::cout << "sent server bytes: " << server_bytes << std::endl;
        return RET_NO_ANS;
    } else {
        size_t num_inputs, num_valid = 0;;
        recv_size(serverfd, num_inputs);
        uint64_t a[B+1];
        memset(a, 0, sizeof(a));
        bool* const valid = new bool[num_inputs];

        for (unsigned int i =0; i < num_inputs; i++) {
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
        recv_in(serverfd, &b[0], share_sz);

        std::cout << "Final valid count: " << num_valid << " / " << total_inputs << std::endl;
        std::cout << "compute time: " << sec_from(start) << std::endl;
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
returnType var_op(const initMsg msg, const int clientfd, const int serverfd, const int server_num, double& ans) {
    auto start = clock_start();

    typedef std::tuple <uint64_t, uint64_t, ClientPacket*> sharetype;
    std::unordered_map<std::string, sharetype> share_map;

    VarShare share;
    const uint64_t max_val = 1ULL << num_bits;
    const unsigned int total_inputs = msg.num_of_inputs;
    const size_t nbits[2] = {num_bits, num_bits * 2};

    // Just for getting sizes
    Circuit* const mock_circuit = CheckVar();
    const size_t NMul = mock_circuit->NumMulGates();
    delete mock_circuit;

    int num_bytes = 0;
    for (unsigned int i = 0; i < total_inputs; i++) {
        num_bytes += recv_in(clientfd, &share, sizeof(VarShare));
        const std::string pk(share.pk, share.pk + PK_LENGTH);

        ClientPacket* packet = new ClientPacket(NMul);
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

        uint64_t* const shares = new uint64_t[2 * num_inputs];
        ClientPacket** const packet = new ClientPacket*[num_inputs];

        std::string* const pk_list = new std::string[num_inputs];
        Circuit** const circuit = new Circuit*[num_inputs];

        size_t idx = 0;
        for (const auto& share : share_map) {
            server_bytes += send_out(serverfd, &share.first[0], PK_LENGTH);
            pk_list[idx] = share.first;
            idx++;
        }
        std::cout << "PK time: " << sec_from(start2) << std::endl;

        for (unsigned int i = 0; i < num_inputs; i++) {
            uint64_t val = 0, val2 = 0;
            std::tie(val, val2, packet[i]) = share_map[pk_list[i]];
            circuit[i] = CheckVar();
            shares[2 * i] = val;
            shares[2 * i + 1] = val2;
        }
        fmpz_t* const shares_p = share_convert(num_inputs, 2,
                                               nbits, shares);
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
        // std::cout << "snip time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        // Convert
        fmpz_t* b = accumulate(num_inputs, 2, shares_p, valid);

        std::cout << "compute time: " << sec_from(start) << std::endl;

        send_fmpz(serverfd, b[0]);
        send_fmpz(serverfd, b[1]);

        delete[] valid;
        clear_fmpz_array(b, 2);
        clear_fmpz_array(shares_p, num_inputs * 2);

        std::cout << "sent non-snip server bytes: " << server_bytes << std::endl;
        return RET_NO_ANS;
    } else {
        size_t num_inputs, num_valid = 0;
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
        }
        std::cout << "PK time: " << sec_from(start2) << std::endl;
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
        const bool* const snip_valid = validate_snips(
            num_inputs, 2, serverfd, server_num, circuit, packet, shares_p);

        for (unsigned int i = 0; i < num_inputs; i++) {
            valid[i] &= snip_valid[i];
            if (valid[i])
                num_valid++;

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
        std::cout << "snip time: " << sec_from(start2) << std::endl;

        // Convert
        fmpz_t* a = accumulate(num_inputs, 2, shares_p, valid);

        std::cout << "compute time: " << sec_from(start) << std::endl;

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

        fmpz_add(b, b, a[0]); fmpz_mod(b, b, Int_Modulus);
        fmpz_add(b2, b2, a[1]); fmpz_mod(b2, b2, Int_Modulus);

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

    // [x], y, [x2], [xy]
    typedef std::tuple <uint64_t*, uint64_t, uint64_t*, uint64_t*, ClientPacket*> sharetype;
    std::unordered_map<std::string, sharetype> share_map;

    const uint64_t max_val = 1ULL << num_bits;
    const unsigned int total_inputs = msg.num_of_inputs;
    size_t nbits[num_fields];
    for (unsigned int i = 0; i < num_fields; i++)
        nbits[i] = num_bits * (i >= degree ? 2 : 1);

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

        ClientPacket* packet = new ClientPacket(NMul);
        int packet_bytes = recv_ClientPacket(clientfd, packet, NMul);
        num_bytes += packet_bytes;

        if ((share_map.find(pk) != share_map.end())
            or (not sizes_valid)
            or (packet_bytes  <= 0)
            ) {
            delete packet;
            continue;
        }

        share_map[pk] = {share.x_vals, share.y, share.x2_vals, share.xy_vals, packet};
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

        uint64_t* const shares = new uint64_t[num_inputs * num_fields];
        ClientPacket** const packet = new ClientPacket*[num_inputs];

        std::string* const pk_list = new std::string[num_inputs];
        Circuit** const circuit = new Circuit*[num_inputs];

        size_t idx = 0;
        for (const auto& share : share_map) {
            server_bytes += send_out(serverfd, &share.first[0], PK_LENGTH);
            pk_list[idx] = share.first;
            idx++;
        }
        std::cout << "PK time: " << sec_from(start2) << std::endl;
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
        std::cout << "snip time: " << sec_from(start2) << std::endl;

        // Convert
        fmpz_t* b = accumulate(num_inputs, num_fields, shares_p, valid);

        std::cout << "compute time: " << sec_from(start) << std::endl;

        for (unsigned int j = 0; j < num_fields; j++)
            send_fmpz(serverfd, b[j]);

        delete[] valid;
        clear_fmpz_array(b, num_fields);
        clear_fmpz_array(shares_p, num_inputs * num_fields);

        std::cout << "sent non-snip server bytes: " << server_bytes << std::endl;

        return RET_NO_ANS;
    } else {
        size_t num_inputs, num_valid = 0;
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
        }
        std::cout << "PK time: " << sec_from(start2) << std::endl;
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
        const bool* const snip_valid = validate_snips(
            num_inputs, num_fields, serverfd, server_num, circuit,
            packet, shares_p);

        for (unsigned int i = 0; i < num_inputs; i++) {
            valid[i] &= snip_valid[i];
            if (valid[i])
                num_valid++;

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
        std::cout << "snip time: " << sec_from(start2) << std::endl;
        start2 = clock_start();

        // Convert
        fmpz_t* a = accumulate(num_inputs, num_fields, shares_p, valid);

        std::cout << "compute time: " << sec_from(start) << std::endl;

        delete[] valid;
        clear_fmpz_array(shares_p, num_inputs * num_fields);


        uint64_t* const x_accum = new uint64_t[degree + num_quad];
        memset(x_accum, 0, (degree + num_quad) * sizeof(uint64_t));
        uint64_t* const y_accum = new uint64_t[degree];
        memset(y_accum, 0, (degree) * sizeof(uint64_t));

        fmpz_t b; fmpz_init(b);

        x_accum[0] = num_valid;
        for (unsigned int j = 0; j < num_x; j++) {
            recv_fmpz(serverfd, b);
            fmpz_add(b, b, a[j]);
            fmpz_mod(b, b, Int_Modulus);
            x_accum[1 + j] = fmpz_get_si(b);
        }
        recv_fmpz(serverfd, b);
        fmpz_add(b, b, a[num_x]);
        fmpz_mod(b, b, Int_Modulus);
        y_accum[0] = fmpz_get_si(b);

        for (unsigned int j = 0; j < num_quad; j++) {
            recv_fmpz(serverfd, b);
            fmpz_add(b, b, a[degree + j]);
            fmpz_mod(b, b, Int_Modulus);
            x_accum[1 + num_x + j] = fmpz_get_si(b);
        }
        for (unsigned int j = 0; j < num_x; j++) {
            recv_fmpz(serverfd, b);
            fmpz_add(b, b, a[degree + num_quad + j]);
            fmpz_mod(b, b, Int_Modulus);
            y_accum[1 + j] = fmpz_get_si(b);
        }

        clear_fmpz_array(a, num_fields);
        fmpz_clear(b);


        std::cout << "Final valid count: " << num_valid << " / " << total_inputs << std::endl;
        std::cout << "sent non-snip server bytes: " << server_bytes << std::endl;
        if (num_valid < total_inputs * (1 - INVALID_THRESHOLD)) {
            std::cout << "Failing, This is less than the invalid threshold of " << INVALID_THRESHOLD << std::endl;
            return RET_INVALID;
        }

        double* c = SolveLinReg(degree, x_accum, y_accum);
        std::cout << "Estimate: y = ";
        for (unsigned int i = 0; i < degree; i++) {
            if (i > 0) std::cout << " + ";
            std::cout << c[i];
            if (i > 0) std::cout << " * x_" << (i-1);
        }
        std::cout << std::endl;
        delete[] x_accum;
        delete[] y_accum;
        delete[] c;

        return RET_ANS;
    }
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cout << "Usage: ./bin/server server_num(0/1) this_client_port server0_port num_bits" << endl;
        return 1;
    }

    const int server_num = atoi(argv[1]);  // Server # 1 or # 2
    const int client_port = atoi(argv[2]); // port of this server, for the client
    const int server_port = atoi(argv[3]); // port of this server, for the other server

    std::cout << "This server is server # " << server_num << std::endl;
    std::cout << "  Listening for client on " << client_port << std::endl;
    std::cout << "  Listening for server on " << server_port << std::endl;

    if (argc >= 5)
        num_bits = atoi(argv[4]);

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

    correlated_store = new CorrelatedStore(serverfd, server_num, ot0, ot1, num_bits, CACHE_SIZE, LAZY_PRECOMPUTE, true, OVER_PRECOMPUTE);

    int sockfd, newsockfd;
    sockaddr_in addr;

    bind_and_listen(addr, sockfd, client_port, 1);

    while(1) {
        // Refresh randomX if used too much
        if (randx_uses > RANDOMX_THRESHOLD) {
            randx_uses = 0;
            sync_randomX(serverfd, server_num, randomX);
            // Update precomps
            for (const auto& pair : precomp_store)
                pair.second -> setCheckerPrecomp(randomX);
        }

        correlated_store->maybeUpdate();

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
            if (ret == RET_ANS) {
                ;
            }

            std::cout << "Total time  : " << sec_from(start) << std::endl;
        } else if (msg.type == NONE_OP) {
            std::cout << "Empty client message" << std::endl;
        } else {
            std::cout << "Unrecognized message type: " << msg.type << std::endl;
        }
        close(newsockfd);
    }

    delete correlated_store;
    for (const auto& precomp : precomp_store)
        delete precomp.second;

    delete ot0;
    delete ot1;
    return 0;
}
