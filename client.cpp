/*
Simulates a group of num_submission clients that communicate with the servers.

General layout:

X_op_helper: Makes then sends a batch of client requests
X_op: Sends init msg, then sends either one batch or a bunch in serial
x_op_invalid: For testing/debugging, does a basic run with intentionally invalid clients.

TODO: invalid linreg
TODO: linreg degree as arg
TODO: probably leave invalid and batch as constants
*/

// TODO: Eventually htonl/ntohl wrappers on shares. Fine when client/server on same machine.

#include "client.h"

#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <unistd.h>

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#include "circuit.h"
#include "net_share.h"
#include "proto.h"
#include "types.h"

// #define SERVER0_IP "52.87.230.64"
// #define SERVER1_IP "54.213.189.18"

#define SERVER0_IP "127.0.0.1"
#define SERVER1_IP "127.0.0.1"

// Whether to include invalid client submissions for testing
#define DEBUG_INVALID false
// Whether to have the client batch or not
#define CLIENT_BATCH true

uint32_t num_bits;
uint64_t max_int;
uint32_t linreg_degree = 2;

int sockfd0, sockfd1;

void error_exit(const char* const msg) {
    perror(msg);
    exit(EXIT_FAILURE);
}

std::string pub_key_to_hex(const uint64_t* const key) {
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(16) << std::hex << key[0];
    ss << std::setfill('0') << std::setw(16) << std::hex << key[1];
    return ss.str();
}

std::string make_pk(emp::PRG prg) {
    emp::block b;
    prg.random_block(&b, 1);
    return pub_key_to_hex((uint64_t*)&b);
}

int send_maxshare(const int server_num, const MaxShare& maxshare, const unsigned int B) {
    const int sock = (server_num == 0) ? sockfd0 : sockfd1;

    int ret = send(sock, (void *)&maxshare, PK_LENGTH, 0);

    for (unsigned int i = 0; i <= B; i++)
        ret += send(sock, (void *)&(maxshare.arr[i]), sizeof(uint32_t), 0);
    return ret;
}

int send_linregshare(const int server_num, const LinRegShare& share,  const size_t degree) {
    const int sock = (server_num == 0) ? sockfd0 : sockfd1;

    const size_t num_x = degree - 1;
    const size_t num_quad = num_x * (num_x + 1) / 2;

    int ret = send(sock, (void*) &share, PK_LENGTH, 0);

    for (unsigned int i = 0; i < num_x; i++)
        ret += send_uint64(sock, share.x_vals[i]);

    ret += send_uint64(sock, share.y);

    for (unsigned int i = 0; i < num_quad; i++)
        ret += send_uint64(sock, share.x2_vals[i]);

    for (unsigned int i = 0; i < num_x; i++)
        ret += send_uint64(sock, share.xy_vals[i]);

    return ret;
}

// Wrapper around send, with error catching.
int send_to_server(const int server, const void* const buffer, const size_t n, const int flags = 0) {
    const int socket = (server == 0 ? sockfd0 : sockfd1);
    int ret = send(socket, buffer, n, flags);
    if (ret < 0) error_exit("Failed to send to server ");
    return ret;
}

int bit_sum_helper(const std::string protocol, const size_t numreqs,
                   unsigned int &ans, const initMsg* const msg_ptr = nullptr) {
    auto start = clock_start();
    int num_bytes = 0;

    bool real_val, share0, share1;

    // Can't use a fixed key, or serial will have the same key every time
    emp::PRG prg;

    BitShare* const bitshare0 = new BitShare[numreqs];
    BitShare* const bitshare1 = new BitShare[numreqs];
    for (unsigned int i = 0; i < numreqs; i++) {
        prg.random_bool(&real_val, 1);
        prg.random_bool(&share0, 1);
        share1 = share0 ^ real_val;
        ans += real_val;

        const std::string pk_s = make_pk(prg);
        const char* const pk = pk_s.c_str();

        // std::cout << pk << ": " << std::noboolalpha << real_vals[i] << " = " << shares0[i] << " ^ " << shares1[i] << std::endl;

        memcpy(bitshare0[i].pk, &pk[0], PK_LENGTH);
        bitshare0[i].val = share0;

        memcpy(bitshare1[i].pk, &pk[0], PK_LENGTH);
        bitshare1[i].val = share1;
    }
    if (numreqs > 1)
        std::cout << "batch make:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;

    start = clock_start();
    if (msg_ptr != nullptr) {
        num_bytes += send_to_server(0, msg_ptr, sizeof(initMsg));
        num_bytes += send_to_server(1, msg_ptr, sizeof(initMsg));
    }
    for (unsigned int i = 0; i < numreqs; i++) {
        num_bytes += send_to_server(0, &bitshare0[i], sizeof(BitShare));
        num_bytes += send_to_server(1, &bitshare1[i], sizeof(BitShare));
    }

    delete[] bitshare0;
    delete[] bitshare1;

    if (numreqs > 1)
        std::cout << "batch send:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;

    return num_bytes;
}

void bit_sum(const std::string protocol, const size_t numreqs) {
    unsigned int ans = 0;
    int num_bytes = 0;
    initMsg msg;
    msg.num_of_inputs = numreqs;
    msg.type = BIT_SUM;

    if (CLIENT_BATCH) {
        num_bytes += bit_sum_helper(protocol, numreqs, ans, &msg);
    } else {
        auto start = clock_start();
        for (unsigned int i = 0; i < numreqs; i++)
            num_bytes += bit_sum_helper(protocol, 1, ans, i == 0 ? &msg : nullptr);
        std::cout << "make+send:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;
    }

    std::cout << "Ans : " << ans << std::endl;
    std::cout << "Total sent bytes: " << num_bytes << std::endl;
}

/* 0: pk mismatch
   2: share0 has same pk as 1
   4: share1 has same pk as 3
*/
void bit_sum_invalid(const std::string protocol, const size_t numreqs) {
    initMsg msg;
    msg.num_of_inputs = numreqs;
    msg.type = BIT_SUM;
    send_to_server(0, &msg, sizeof(initMsg));
    send_to_server(1, &msg, sizeof(initMsg));

    emp::block* const b = new block[numreqs];
    bool real_vals[numreqs];
    bool shares0[numreqs];
    bool shares1[numreqs];

    emp::PRG prg(fix_key);
    prg.random_block(b, numreqs);
    prg.random_bool(real_vals, numreqs);
    prg.random_bool(shares0, numreqs);

    int ans = 0;
    std::string pk_str = "";
    for (unsigned int i = 0; i < numreqs; i++) {
        BitShare share0, share1;
        const char* prev_pk = pk_str.c_str();
        pk_str = pub_key_to_hex((uint64_t*)&b[i]);
        const char* const pk = pk_str.c_str();

        shares1[i] = real_vals[i]^shares0[i];
        std::cout << i << ": " << std::boolalpha << shares0[i] << " ^ " << shares1[i] << " = " << real_vals[i];
        if (i == 0 or i == 2 or i == 4) {
            std::cout << " (invalid)" << std::endl;
        } else {
            std::cout << std::endl;
            ans += (real_vals[i] ? 1 : 0);
        }

        memcpy(share0.pk, &pk[0], PK_LENGTH);
        share0.val = shares0[i];
        if (i == 0)
            share0.pk[0] = 'q';
        if (i == 2)
            memcpy(share0.pk, &prev_pk[0], PK_LENGTH);

        memcpy(share1.pk, &pk[0], PK_LENGTH);
        share1.val = shares1[i];
        if (i == 4)
            memcpy(share1.pk, &prev_pk[0], PK_LENGTH);

        send_to_server(0, &share0, sizeof(BitShare));
        send_to_server(1, &share1, sizeof(BitShare));
    }
    std::cout << "Ans : " << ans << std::endl;

    delete[] b;
}

int int_sum_helper(const std::string protocol, const size_t numreqs,
                   uint64_t &ans, const initMsg* const msg_ptr = nullptr) {
    auto start = clock_start();
    int num_bytes = 0;

    uint64_t real_val, share0, share1;

    emp::PRG prg;

    IntShare* const intshare0 = new IntShare[numreqs];
    IntShare* const intshare1 = new IntShare[numreqs];

    for (unsigned int i = 0; i < numreqs; i++) {
        prg.random_data(&real_val, sizeof(uint64_t));
        prg.random_data(&share0, sizeof(uint64_t));
        real_val = real_val % max_int;
        share0 = share0 % max_int;
        share1 = share0 ^ real_val;
        ans += real_val;

        const std::string pk_s = make_pk(prg);
        const char* const pk = pk_s.c_str();

        memcpy(intshare0[i].pk, &pk[0], PK_LENGTH);
        intshare0[i].val = share0;

        memcpy(intshare1[i].pk, &pk[0], PK_LENGTH);
        intshare1[i].val = share1;
    }
    if (numreqs > 1)
        std::cout << "batch make:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;

    start = clock_start();
    if (msg_ptr != nullptr) {
        num_bytes += send_to_server(0, msg_ptr, sizeof(initMsg));
        num_bytes += send_to_server(1, msg_ptr, sizeof(initMsg));
    }
    for (unsigned int i = 0; i < numreqs; i++) {
        num_bytes += send_to_server(0, &intshare0[i], sizeof(IntShare));
        num_bytes += send_to_server(1, &intshare1[i], sizeof(IntShare));
    }
    delete[] intshare0;
    delete[] intshare1;

    if (numreqs > 1)
        std::cout << "batch send:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;

    return num_bytes;
}

void int_sum(const std::string protocol, const size_t numreqs) {
    uint64_t ans = 0;
    int num_bytes = 0;
    initMsg msg;
    msg.num_of_inputs = numreqs;
    msg.type = INT_SUM;

    if (CLIENT_BATCH) {
        num_bytes += int_sum_helper(protocol, numreqs, ans, &msg);
    } else {
        auto start = clock_start();
        for (unsigned int i = 0; i < numreqs; i++)
            num_bytes += int_sum_helper(protocol, 1, ans, i == 0 ? &msg : nullptr);
        std::cout << "make+send:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;
    }

    std::cout << "Ans : " << ans << std::endl;
    std::cout << "Total sent bytes: " << num_bytes << std::endl;
}

/* 0: x > max
   1: x share > max
   2: pk mismatch
   4: share0 has same pk as 3
   6: share1 has same pk as 5
*/
void int_sum_invalid(const std::string protocol, const size_t numreqs) {
    initMsg msg;
    msg.num_of_inputs = numreqs;
    msg.type = INT_SUM;
    send_to_server(0, &msg, sizeof(initMsg));
    send_to_server(1, &msg, sizeof(initMsg));

    emp::block* const b = new block[numreqs];
    uint64_t real_vals[numreqs];
    uint64_t shares0[numreqs];
    uint64_t shares1[numreqs];
    uint64_t ans = 0;

    emp::PRG prg(fix_key);
    prg.random_block(b, numreqs);
    prg.random_data(real_vals, numreqs * sizeof(uint64_t));
    prg.random_data(shares0, numreqs * sizeof(uint64_t));

    std::string pk_str = "";
    for (unsigned int i = 0; i < numreqs; i++) {
        if (i != 0)
            real_vals[i] = real_vals[i] % max_int;
        if (i != 1)
            shares0[i] = shares0[i] % max_int;
        shares1[i] = real_vals[i] ^ shares0[i];
        std::cout << "real_vals[" << i << "] = " << real_vals[i] << " = " << shares0[i] << " ^ " << shares1[i];
        if (i <= 2 or i == 4 or i == 6) {
            std::cout << " (invalid)" << std::endl;
        } else {
            std::cout << std::endl;
            ans += real_vals[i];
        }

        IntShare share0, share1;
        const char* prev_pk = pk_str.c_str();
        pk_str = pub_key_to_hex((uint64_t*)&b[i]);
        const char* const pk = pk_str.c_str();

        memcpy(share0.pk, &pk[0], PK_LENGTH);
        share0.val = shares0[i];
        if (i == 2)
            share0.pk[0] = 'q';
        if (i == 4)
            memcpy(share0.pk, &prev_pk[0], PK_LENGTH);

        memcpy(share1.pk, &pk[0], PK_LENGTH);
        share1.val = shares1[i];
        if (i == 6)
            memcpy(share1.pk, &prev_pk[0], PK_LENGTH);

        send_to_server(0, &share0, sizeof(IntShare));
        send_to_server(1, &share1, sizeof(IntShare));
    }

    std::cout << "Ans : " << ans << std::endl;

    delete[] b;
}

int xor_op_helper(const std::string protocol, const size_t numreqs,
                  bool &ans, const initMsg* const msg_ptr = nullptr) {
    auto start = clock_start();
    int num_bytes = 0;

    bool value;
    uint64_t encoded, share0, share1;

    emp::PRG prg;

    IntShare* const intshare0 = new IntShare[numreqs];
    IntShare* const intshare1 = new IntShare[numreqs];

    for (unsigned int i = 0; i < numreqs; i++) {
        prg.random_bool(&value, 1);
        if (protocol == "ANDOP") {
            ans &= value;
            if (value)
                encoded = 0;
            else
                prg.random_data(&encoded, sizeof(uint64_t));
        } else if (protocol == "OROP") {
            ans |= value;
            if (not value)
                encoded = 0;
            else
                prg.random_data(&encoded, sizeof(uint64_t));
        }
        prg.random_data(&share0, sizeof(uint64_t));
        share1 = share0 ^ encoded;
    
        const std::string pk_s = make_pk(prg);
        const char* const pk = pk_s.c_str();

        memcpy(intshare0[i].pk, &pk[0], PK_LENGTH);
        intshare0[i].val = share0;

        memcpy(intshare1[i].pk, &pk[0], PK_LENGTH);
        intshare1[i].val = share1;
    }
    if (numreqs > 1)
        std::cout << "batch make:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;
    start = clock_start();
    if (msg_ptr != nullptr) {
        num_bytes += send_to_server(0, msg_ptr, sizeof(initMsg));
        num_bytes += send_to_server(1, msg_ptr, sizeof(initMsg));
    }
    for (unsigned int i = 0; i < numreqs; i++) {
        num_bytes += send_to_server(0, &intshare0[i], sizeof(IntShare));
        num_bytes += send_to_server(1, &intshare1[i], sizeof(IntShare));
    }

    delete[] intshare0;
    delete[] intshare1;

    if (numreqs > 1)
        std::cout << "batch send:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;

    return num_bytes;
}

void xor_op(const std::string protocol, const size_t numreqs) {
    bool ans;
    int num_bytes = 0;
    initMsg msg;
    msg.num_of_inputs = numreqs;
    if (protocol == "ANDOP") {
        msg.type = AND_OP;
        ans = true;
    } else if (protocol == "OROP") {
        msg.type = OR_OP;
        ans = false;
    } else {
        return;
    }

    if (CLIENT_BATCH) {
        num_bytes += xor_op_helper(protocol, numreqs, ans, &msg);
    } else {
        auto start = clock_start();
        for (unsigned int i = 0; i < numreqs; i++)
            num_bytes += xor_op_helper(protocol, 1, ans, i == 0 ? &msg : nullptr);
        std::cout << "make+send:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;
    }

    std::cout << "Ans : " << std::boolalpha << ans << std::endl;
    std::cout << "Total sent bytes: " << num_bytes << std::endl;
}

/* 0: pk mismatch
   2: share0 has same pk as 1
   4: share1 has same pk as 3
*/
void xor_op_invalid(const std::string protocol, const size_t numreqs) {
    initMsg msg;
    msg.num_of_inputs = numreqs;
    bool ans;
    if (protocol == "ANDOP") {
        msg.type = AND_OP;
        ans = true;
    } else if (protocol == "OROP") {
        msg.type = OR_OP;
        ans = false;
    } else {
        return;
    }
    send_to_server(0, &msg, sizeof(initMsg));
    send_to_server(1, &msg, sizeof(initMsg));

    emp::block* const b = new block[numreqs];
    bool values[numreqs];
    uint64_t encoded_values[numreqs];
    uint64_t shares0[numreqs];
    uint64_t shares1[numreqs];

    emp::PRG prg(fix_key);

    prg.random_block(b, numreqs);
    prg.random_bool(values, numreqs);
    prg.random_data(encoded_values, numreqs * sizeof(uint64_t));
    prg.random_data(shares0, numreqs * sizeof(uint64_t));

    for (unsigned int i = 0; i < numreqs; i++) {
        std::cout << "val[" << i << "] = " << std::boolalpha << values[i];
        if (i == 0 or i == 2 or i == 4) {
            std::cout << " (invalid) " << std::endl;
            continue;
        }
        std::cout << std::endl;
        if (protocol == "ANDOP")
            ans &= values[i];
        if (protocol == "OROP")
            ans |= values[i];
    }

    // encode step. set to all 0's for values that don't force the ans.
    if (protocol == "ANDOP")
        for (unsigned int i = 0; i < numreqs; i++)
            if (values[i])
                encoded_values[i] = 0;
    if (protocol == "OROP")
        for (unsigned int i = 0; i < numreqs; i++)
            if (!values[i])
                encoded_values[i] = 0;

    std::string pk_str = "";

    // Share splitting. Same as int sum. Sum of shares = encoded value
    for (unsigned int i = 0; i < numreqs; i++) {
        shares1[i] = encoded_values[i] ^ shares0[i];

        IntShare share0, share1;
        const char* prev_pk = pk_str.c_str();
        pk_str = pub_key_to_hex((uint64_t*)&b[i]);
        const char* const pk = pk_str.c_str();

        memcpy(share0.pk, &pk[0], PK_LENGTH);
        share0.val = shares0[i];
        if (i == 0)
            share0.pk[0] = 'q';
        if (i == 2)
            memcpy(share0.pk, &prev_pk[0], PK_LENGTH);

        memcpy(share1.pk, &pk[0], PK_LENGTH);
        share1.val = shares1[i];
        if (i == 4)
            memcpy(share1.pk, &prev_pk[0], PK_LENGTH);

        send_to_server(0, &share0, sizeof(IntShare));
        send_to_server(1, &share1, sizeof(IntShare));
    }

    std::cout << "Ans : " << std::boolalpha << ans << std::endl;

    delete[] b;
}

int max_op_helper(const std::string protocol, const size_t numreqs,
                  const unsigned int B, uint32_t &ans,
                  const initMsg* const msg_ptr = nullptr) {
    auto start = clock_start();
    int num_bytes = 0;

    start = clock_start();

    uint32_t value;
    uint32_t* const or_encoded_array = new uint32_t[B+1];
    uint32_t* const share0 = new uint32_t[B+1];
    uint32_t* const share1 = new uint32_t[B+1];

    emp::PRG prg;

    MaxShare* const maxshare0 = new MaxShare[numreqs];
    MaxShare* const maxshare1 = new MaxShare[numreqs];
    for (unsigned int i = 0; i < numreqs; i++) {
        prg.random_data(&value, sizeof(uint32_t));
        value = value % (B + 1);

        if (protocol == "MAXOP")
            ans = (value > ans ? value : ans);
        if (protocol == "MINOP")
            ans = (value < ans ? value : ans);

        prg.random_data(or_encoded_array, (B+1)*sizeof(uint32_t));
        prg.random_data(share0, (B+1)*sizeof(uint32_t));

        uint32_t v = 0;
        if (protocol == "MAXOP")
            v = value;
        if (protocol == "MINOP")
            v = B - value;

        for (unsigned int j = v + 1; j <= B ; j++)
            or_encoded_array[j] = 0;

        for (unsigned int j = 0; j <= B; j++)
            share1[j] = share0[j] ^ or_encoded_array[j];

        const std::string pk_s = make_pk(prg);
        const char* const pk = pk_s.c_str();

        memcpy(maxshare0[i].pk, &pk[0], PK_LENGTH);
        maxshare0[i].arr = new uint32_t[B+1];
        memcpy(maxshare0[i].arr, share0, (B+1)*sizeof(uint32_t));

        memcpy(maxshare1[i].pk, &pk[0], PK_LENGTH);
        maxshare1[i].arr = new uint32_t[B+1];
        memcpy(maxshare1[i].arr, share1, (B+1)*sizeof(uint32_t));
    }
    delete[] or_encoded_array;
    delete[] share0;
    delete[] share1;
    if (numreqs > 1)
        std::cout << "batch make:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;

    start = clock_start();
    if (msg_ptr != nullptr) {
        num_bytes += send_to_server(0, msg_ptr, sizeof(initMsg));
        num_bytes += send_to_server(1, msg_ptr, sizeof(initMsg));
    }
    for (unsigned int i = 0; i < numreqs; i++) {
        num_bytes += send_maxshare(0, maxshare0[i], B);
        num_bytes += send_maxshare(1, maxshare1[i], B);

        delete[] maxshare0[i].arr;
        delete[] maxshare1[i].arr;
    }

    delete[] maxshare0;
    delete[] maxshare1;

    if (numreqs > 1)
        std::cout << "batch send:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;

    return num_bytes;
}

void max_op(const std::string protocol, const size_t numreqs) {
    const uint64_t B = max_int;

    uint32_t ans;
    int num_bytes = 0;
    initMsg msg;
    msg.num_of_inputs = numreqs;
    msg.max_inp = B;
    if (protocol == "MAXOP") {
        msg.type = MAX_OP;
        ans = 0;
    } else if (protocol == "MINOP") {
        msg.type = MIN_OP;
        ans = B;
    } else {
        return;
    }

    if (CLIENT_BATCH) {
        num_bytes += max_op_helper(protocol, numreqs, B, ans, &msg);
    } else {
        auto start = clock_start();
        for (unsigned int i = 0; i < numreqs; i++)
            num_bytes += max_op_helper(protocol, 1, B, ans, i == 0 ? &msg : nullptr);
        std::cout << "make+send:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;
    }

    std::cout << "Ans : " << ans << std::endl;
    std::cout << "Total sent bytes: " << num_bytes << std::endl;
}

/* 0: pk mismatch
   2: share0 has same pk as 1
   4: share1 has same pk as 3
*/
void max_op_invalid(const std::string protocol, const size_t numreqs) {
    const unsigned int B = 250;
    initMsg msg;
    msg.num_of_inputs = numreqs;
    msg.max_inp = B;
    emp::PRG prg(fix_key);
    uint32_t ans;
    if (protocol == "MAXOP") {
        msg.type = MAX_OP;
        ans = 0;
    } else if (protocol == "MINOP") {
        msg.type = MIN_OP;
        ans = B;
    } else {
        return;
    }

    emp::block* const b = new block[numreqs];
    prg.random_block(b, numreqs);

    uint32_t values[numreqs];
    uint32_t or_encoded_array[B+1];
    uint32_t shares0[B+1];
    uint32_t shares1[B+1];
    prg.random_data(values, numreqs * sizeof(uint32_t));

    send_to_server(0, &msg, sizeof(initMsg), 0);
    send_to_server(1, &msg, sizeof(initMsg), 0);

    std::string pk_str = "";

    for (unsigned int i = 0; i < numreqs; i++) {
        MaxShare share0, share1;
        values[i] = values[i] % (B + 1);
        std::cout << "value[" << i << "] = " << values[i];
        if (i == 0 or i == 2 or i == 4) {
            std::cout << " (invalid)" << std::endl;
        } else {
            std::cout << std::endl;
            if (protocol == "MAXOP")
                ans = (values[i] > ans? values[i] : ans);
            if (protocol == "MINOP")
                ans = (values[i] < ans? values[i] : ans);
        }

        prg.random_data(or_encoded_array, (B+1)*sizeof(uint32_t));
        prg.random_data(shares0, (B+1)*sizeof(uint32_t));

        // min(x) = -max(-x) = B - max(B - x)
        uint32_t v = 0;
        if (protocol == "MAXOP")
            v = values[i];
        if (protocol == "MINOP")
            v = B - values[i];

        for (unsigned int j = v + 1; j <= B ; j++)
            or_encoded_array[j] = 0;

        for (unsigned int j = 0; j <= B; j++)
            shares1[j] = shares0[j] ^ or_encoded_array[j];

        const char* prev_pk = pk_str.c_str();
        pk_str = pub_key_to_hex((uint64_t*)&b[i]);
        const char* const pk = pk_str.c_str();

        memcpy(share0.pk, &pk[0], PK_LENGTH);
        share0.arr = shares0;
        if (i == 0)
            share0.pk[0] = 'q';
        if (i == 2)
            memcpy(share0.pk, &prev_pk[0], PK_LENGTH);

        memcpy(share1.pk, &pk[0], PK_LENGTH);
        share1.arr = shares1;
        if (i == 4)
            memcpy(share1.pk, &prev_pk[0], PK_LENGTH);

        send_maxshare(0, share0, B);
        send_maxshare(1, share1, B);
    }

    std::cout << "Ans : " << ans << std::endl;

    delete[] b;
}

int var_op_helper(const std::string protocol, const size_t numreqs,
                  uint64_t& sum, uint64_t& sumsquared,
                  const initMsg* const msg_ptr = nullptr) {
    auto start = clock_start();
    int num_bytes = 0;

    uint64_t real_val, share0, share1, share0_2, share1_2;

    fmpz_t inp[2];
    fmpz_init(inp[0]);
    fmpz_init(inp[1]);

    emp::PRG prg;

    VarShare* const varshare0 = new VarShare[numreqs];
    VarShare* const varshare1 = new VarShare[numreqs];
    ClientPacket** const packet0 = new ClientPacket*[numreqs];
    ClientPacket** const packet1 = new ClientPacket*[numreqs];

    for (unsigned int i = 0; i < numreqs; i++) {
        prg.random_data(&real_val, sizeof(uint64_t));
        prg.random_data(&share0, sizeof(uint64_t));
        prg.random_data(&share0_2, sizeof(uint64_t));
        real_val = real_val % max_int;
        share0 = share0 % max_int;
        share1 = share0 ^ real_val;
        const uint64_t squared = real_val * real_val;
        share0_2 = share0_2 % (max_int * max_int);
        share1_2 = share0_2 ^ squared;
        sum += real_val;
        sumsquared += squared;

        const std::string pk_s = make_pk(prg);
        const char* const pk = pk_s.c_str();

        memcpy(varshare0[i].pk, &pk[0], PK_LENGTH);
        varshare0[i].val = share0;
        varshare0[i].val_squared = share0_2;

        memcpy(varshare1[i].pk, &pk[0], PK_LENGTH);
        varshare1[i].val = share1;
        varshare1[i].val_squared = share1_2;

        fmpz_set_si(inp[0], real_val);
        fmpz_set_si(inp[1], squared);
        Circuit* const circuit = CheckVar();
        circuit->Eval(inp);
        packet0[i] = new ClientPacket(circuit->N(), circuit->NumMulInpGates());
        packet1[i] = new ClientPacket(circuit->N(), circuit->NumMulInpGates());
        share_polynomials(circuit, packet0[i], packet1[i]);
        delete circuit;
    }
    if (numreqs > 1)
        std::cout << "batch make:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;

    start = clock_start();
    if (msg_ptr != nullptr) {
        num_bytes += send_to_server(0, msg_ptr, sizeof(initMsg));
        num_bytes += send_to_server(1, msg_ptr, sizeof(initMsg));
    }
    for (unsigned int i = 0; i < numreqs; i++) {
        num_bytes += send_to_server(0, &varshare0[i], sizeof(VarShare));
        num_bytes += send_to_server(1, &varshare1[i], sizeof(VarShare));

        num_bytes += send_ClientPacket(sockfd0, packet0[i]);
        num_bytes += send_ClientPacket(sockfd1, packet1[i]);

        delete packet0[i];
        delete packet1[i];
    }

    delete[] varshare0;
    delete[] varshare1;
    delete[] packet0;
    delete[] packet1;
    fmpz_clear(inp[0]);
    fmpz_clear(inp[1]);

    if (numreqs > 1)
        std::cout << "batch send:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;

    return num_bytes;
}

void var_op(const std::string protocol, const size_t numreqs) {
    if (num_bits > 31)
        error_exit("Num bits is too large. x^2 > 2^64.");

    uint64_t sum = 0, sumsquared = 0;
    int num_bytes = 0;
    initMsg msg;
    msg.num_of_inputs = numreqs;
    if (protocol == "VAROP") {
        msg.type = VAR_OP;
    } else if (protocol == "STDDEVOP") {
        msg.type = STDDEV_OP;
    } else {
        return;
    }

    if (CLIENT_BATCH) {
        num_bytes += var_op_helper(protocol, numreqs, sum, sumsquared, &msg);
    } else {
        auto start = clock_start();
        for (unsigned int i = 0; i < numreqs; i++)
            num_bytes += var_op_helper(protocol, 1, sum, sumsquared, i == 0 ? &msg : nullptr);
        std::cout << "make+send:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;
    }

    const double ex = 1. * sum / numreqs;
    const double ex2 = 1. * sumsquared / numreqs;
    double ans = ex2 - (ex * ex);
    std::cout << "E[X^2] - E[X]^2 = " << ex2 << " - (" << ex << ")^2 = " << ans << std::endl;
    if (protocol == "STDDEVOP")
        ans = sqrt(ans);
    std::cout << "True Ans: " << ans << std::endl;
    std::cout << "Total sent bytes: " << num_bytes << std::endl;
}

/* 0: x > max
   1: x shares > max
   2: x^2 shares > max
   3: x^2 = x * x + junk, run through snip
   4: x^2 = x * x + junk, run snip with correct x^2  NOT CAUGHT, Server doesn't check equality
   5: p0 N corruption
   6: p1 const value corruption
   7: p0 wire corruption
   8: p1 triple corruption
   9: pk mismatch
   11: share0 has same pk as 10
   13: share1 has same pk as 12
*/
void var_op_invalid(const std::string protocol, const size_t numreqs) {
    initMsg msg;
    msg.num_of_inputs = numreqs;
    if (protocol == "VAROP") {
        msg.type = VAR_OP;
    } else if (protocol == "STDDEVOP") {
        msg.type = STDDEV_OP;
    } else {
        return;
    }
    send_to_server(0, &msg, sizeof(initMsg));
    send_to_server(1, &msg, sizeof(initMsg));

    emp::block* const b = new block[numreqs];
    uint64_t real_vals[numreqs];
    // shares of x
    uint64_t shares0[numreqs];
    uint64_t shares1[numreqs];
    // shares of x^2
    uint64_t shares0_squared[numreqs];
    uint64_t shares1_squared[numreqs];
    uint64_t sum = 0, sumsquared = 0, numvalid = 0;

    emp::PRG prg(fix_key);
    prg.random_block(b, numreqs);
    prg.random_data(real_vals, numreqs * sizeof(uint64_t));
    prg.random_data(shares0, numreqs * sizeof(uint64_t));
    prg.random_data(shares0_squared, numreqs * sizeof(uint64_t));

    std::cout << "max_int: " << max_int << std::endl;

    fmpz_t inp[2];
    fmpz_init(inp[0]);
    fmpz_init(inp[1]);

    std::string pk_str = "";

    for (unsigned int i = 0; i < numreqs; i++) {
        real_vals[i] = real_vals[i] % max_int;
        if (i == 0)  // x over cap
            real_vals[i] += max_int;
        shares0[i] = shares0[i] % max_int;
        if (i == 1)  // x shares over capped
            shares0[i] += max_int;
        shares1[i] = real_vals[i] ^ shares0[i];
        uint64_t squared = real_vals[i] * real_vals[i];
        if (i == 3 or i == 4)  // x^2 != x * x
            squared = (squared + 10) % (max_int * max_int);
        if (i != 2)  // x^2 not capped
            shares0_squared[i] = shares0_squared[i] % max_int;
        shares1_squared[i] = squared ^ shares0_squared[i];

        std::cout << i << ": " << real_vals[i] << " = " << shares0[i] << "^" << shares1[i];
        std::cout << ", " << squared << " = " << shares0_squared[i] << "^" << shares1_squared[i];
        if (i <= 9 or i == 11 or i == 13) {
            std::cout << " (invalid)" << std::endl;
        } else {
            std::cout << std::endl;
            sum += real_vals[i];
            sumsquared += squared;
            numvalid++;
        }

        VarShare share0, share1;
        const char* prev_pk = pk_str.c_str();
        pk_str = pub_key_to_hex((uint64_t*)&b[i]);
        const char* const pk = pk_str.c_str();

        memcpy(share0.pk, &pk[0], PK_LENGTH);
        share0.val = shares0[i];
        share0.val_squared = shares0_squared[i];
        if (i == 9)
            share0.pk[0] = 'q';
        if (i == 11)
            memcpy(share0.pk, &prev_pk[0], PK_LENGTH);

        memcpy(share1.pk, &pk[0], PK_LENGTH);
        share1.val = shares1[i];
        share1.val_squared = shares1_squared[i];
        if (i == 13)
            memcpy(share1.pk, &prev_pk[0], PK_LENGTH);

        send_to_server(0, &share0, sizeof(VarShare));
        send_to_server(1, &share1, sizeof(VarShare));
        // SNIP: proof that x^2 = x_squared
        fmpz_set_si(inp[0], real_vals[i]);
        fmpz_set_si(inp[1], real_vals[i] * real_vals[i]);
        if (i == 3)
            fmpz_set_si(inp[1], (real_vals[i] * real_vals[i] + 10) % max_int);
        Circuit* const circuit = CheckVar();
        // Run through circuit to set wires
        circuit->Eval(inp);
        ClientPacket* p0;
        if (i == 5) {
            p0 = new ClientPacket(1, circuit->NumMulInpGates());
        } else {
            p0 = new ClientPacket(circuit->N(), circuit->NumMulInpGates());
        }
        ClientPacket* p1 = new ClientPacket(circuit->N(), circuit->NumMulInpGates());
        share_polynomials(circuit, p0, p1);
        delete circuit;
        if (i == 6)
            fmpz_add_si(p1->f0_s, p1->f0_s, 1);
        if (i == 7)
            fmpz_add_si(p0->WireShares[0], p0->WireShares[0], 1);
        if (i == 8)
            fmpz_add_si(p1->triple_share->shareA, p1->triple_share->shareA, 1);
        send_ClientPacket(sockfd0, p0);
        send_ClientPacket(sockfd1, p1);
        delete p0;
        delete p1;
    }

    std::cout << "sum: " << sum << std::endl;
    std::cout << "sumsquared: " << sumsquared << std::endl;
    std::cout << "numvalid: " << numvalid << std::endl;

    const double ex = 1. * sum / numvalid;
    const double ex2 = 1. * sumsquared / numvalid;
    double ans = ex2 - (ex * ex);
    std::cout << "E[X^2] - E[X]^2 = " << ex2 << " - (" << ex << ")^2 = " << ans << std::endl;
    if (protocol == "STDDEVOP")
        ans = sqrt(ans);
    std::cout << "True Ans: " << ans << std::endl;

    delete[] b;
    fmpz_clear(inp[0]);
    fmpz_clear(inp[1]);
}

// Currently just degree 2
int lin_reg_helper(const std::string protocol, const size_t numreqs,
                   const size_t degree,
                   uint64_t* const x_accum, uint64_t* const y_accum,
                   const initMsg* const msg_ptr = nullptr) {
    auto start = clock_start();

    const size_t num_x = degree - 1;
    const size_t num_quad = num_x * (num_x + 1) / 2;
    const size_t num_fields = 2 * num_x + 1 + num_quad;

    int num_bytes = 0;

    uint64_t* const x_real = new uint64_t[num_x];
    uint64_t* const x_share0 = new uint64_t[num_x];
    uint64_t* const x_share1 = new uint64_t[num_x];
    uint64_t y_real, y_share0, y_share1;
    uint64_t* const x2_real = new uint64_t[num_quad];
    uint64_t* const x2_share0 = new uint64_t[num_quad];
    uint64_t* const x2_share1 = new uint64_t[num_quad];
    uint64_t* const xy_real = new uint64_t[num_x];
    uint64_t* const xy_share0 = new uint64_t[num_x];
    uint64_t* const xy_share1 = new uint64_t[num_x];

    emp::PRG prg;

    LinRegShare* const linshare0 = new LinRegShare[numreqs];
    LinRegShare* const linshare1 = new LinRegShare[numreqs];
    ClientPacket** const packet0 = new ClientPacket*[numreqs];
    ClientPacket** const packet1 = new ClientPacket*[numreqs];

    fmpz_t* inp; new_fmpz_array(&inp, num_fields);

    for (unsigned int i = 0; i < numreqs; i++) {
        prg.random_data(x_real, num_x * sizeof(uint64_t));
        prg.random_data(x_share0, num_x * sizeof(uint64_t));
        prg.random_data(x2_share0, num_quad * sizeof(uint64_t));
        prg.random_data(xy_share0, num_x * sizeof(uint64_t));
        prg.random_data(&y_real, sizeof(uint64_t));
        prg.random_data(&y_share0, sizeof(uint64_t));

        for (unsigned int j = 0; j < num_x; j++) {
            x_real[j] = x_real[j] % max_int;
            x_share0[j] = x_share0[j] % max_int;
            x_share1[j] = x_share0[j] ^ x_real[j];
            x_accum[j + 1] += x_real[j];
            // std::cout << "x[" << j << "]  = " << x_real[j] << " = " << x_share0[j] << " ^ " << x_share1[j] << std::endl;
        }
        y_real = y_real % max_int;
        y_share0 = y_share0 % max_int;
        y_share1 = y_share0 ^ y_real;
        y_accum[0] += y_real;
        // std::cout << "y     = " << y_real << " = " << y_share0 << " ^ " << y_share1 << std::endl;

        // std::cout << "x: ";
        // for (unsigned int j = 0; j < num_x; j++) {
        //     if (j > 0) std::cout << ", ";
        //     std::cout << x_real[j];
        // }
        // std::cout << ", y: " << y_real << std::endl;

        unsigned int idx = 0;
        for (unsigned int j = 0; j < num_x; j++) {
            for (unsigned int k = j; k < num_x; k++) {
                x2_real[idx] = x_real[j] * x_real[k];
                x2_share0[idx] = x2_share0[idx] % (max_int * max_int);
                x2_share1[idx] = x2_share0[idx] ^ x2_real[idx];
                // std::cout << "x2[" << idx << "] = " << x2_real[idx] << " = " << x2_share0[idx] << " ^ " << x2_share1[idx] << std::endl;

                x_accum[idx + num_x + 1] += x2_real[idx];
                idx++;
            }
            xy_real[j] = x_real[j] * y_real;
            xy_share0[j] = xy_share0[j] % (max_int * max_int);
            xy_share1[j] = xy_share0[j] ^ xy_real[j];
            // std::cout << "xy[" << j << "] = " << xy_real[j] << " = " << xy_share0[j] << " ^ " << xy_share1[j] << std::endl;
            
            y_accum[j + 1] += xy_real[j];
        }

        // build shares
        const std::string pk_s = make_pk(prg);
        const char* const pk = pk_s.c_str();

        memcpy(linshare0[i].pk, &pk[0], PK_LENGTH);
        linshare0[i].y = y_share0;
        linshare0[i].x_vals = new uint64_t[num_x];
        linshare0[i].x2_vals = new uint64_t[num_quad];
        linshare0[i].xy_vals = new uint64_t[num_x];
        memcpy(linshare0[i].x_vals, x_share0, num_x * sizeof(uint64_t));
        memcpy(linshare0[i].x2_vals, x2_share0, num_quad * sizeof(uint64_t));
        memcpy(linshare0[i].xy_vals, xy_share0, num_x * sizeof(uint64_t));

        memcpy(linshare1[i].pk, &pk[0], PK_LENGTH);
        linshare1[i].y = y_share1;
        linshare1[i].x_vals = new uint64_t[num_x];
        linshare1[i].x2_vals = new uint64_t[num_quad];
        linshare1[i].xy_vals = new uint64_t[num_x];
        memcpy(linshare1[i].x_vals, x_share1, num_x * sizeof(uint64_t));
        memcpy(linshare1[i].x2_vals, x2_share1, num_quad * sizeof(uint64_t));
        memcpy(linshare1[i].xy_vals, xy_share1, num_x * sizeof(uint64_t));

        // build packets
        for (unsigned int j = 0; j < num_x; j++)
            fmpz_set_si(inp[j], x_real[j]);
        fmpz_set_si(inp[num_x], y_real);
        for (unsigned int j = 0; j < num_quad; j++)
            fmpz_set_si(inp[j + num_x + 1], x2_real[j]);
        for (unsigned int j = 0; j < num_x; j++)
            fmpz_set_si(inp[j + num_x + num_quad + 1], xy_real[j]);

        // for (unsigned int j = 0; j < num_fields; j++) {
        //     std::cout << "inp[" << j << "] = ";
        //     fmpz_print(inp[j]);
        //     std::cout << std::endl;
        // }

        Circuit* const circuit = CheckLinReg(degree);
        circuit->Eval(inp);
        packet0[i] = new ClientPacket(circuit->N(), circuit->NumMulInpGates());
        packet1[i] = new ClientPacket(circuit->N(), circuit->NumMulInpGates());
        share_polynomials(circuit, packet0[i], packet1[i]);
        delete circuit;
    }
    x_accum[0] += numreqs;
    delete[] x_real;
    delete[] x_share0;
    delete[] x_share1;
    delete[] x2_real;
    delete[] x2_share0;
    delete[] x2_share1;
    delete[] xy_real;
    delete[] xy_share0;
    delete[] xy_share1;
    if (numreqs > 1)
        std::cout << "batch make:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;

    start = clock_start();
    if (msg_ptr != nullptr) {
        num_bytes += send_to_server(0, msg_ptr, sizeof(initMsg));
        num_bytes += send_to_server(1, msg_ptr, sizeof(initMsg));
        
        num_bytes += send_size(sockfd0, degree);
        num_bytes += send_size(sockfd1, degree);
    }
    for (unsigned int i = 0; i < numreqs; i++) {
        num_bytes += send_linregshare(0, linshare0[i], degree);
        num_bytes += send_linregshare(1, linshare1[i], degree);

        num_bytes += send_ClientPacket(sockfd0, packet0[i]);
        num_bytes += send_ClientPacket(sockfd1, packet1[i]);

        delete[] linshare0[i].x_vals;
        delete[] linshare0[i].x2_vals;
        delete[] linshare0[i].xy_vals;
        delete[] linshare1[i].x_vals;
        delete[] linshare1[i].x2_vals;
        delete[] linshare1[i].xy_vals;
        delete packet0[i];
        delete packet1[i];
    }

    delete[] linshare0;
    delete[] linshare1;
    delete[] packet0;
    delete[] packet1;
    clear_fmpz_array(inp, num_fields);

    if (numreqs > 1)
        std::cout << "batch send:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;

    return num_bytes;
}

void lin_reg(const std::string protocol, const size_t numreqs) {
    if (num_bits > 31)
        error_exit("Num bits is too large. x^2 > 2^64.");

    const size_t degree = linreg_degree;
    const size_t num_x = degree - 1;
    const size_t num_quad = num_x * (num_x + 1) / 2;
    // const size_t num_fields = 2 * num_x + 1 + num_quad;

    std::cout << "Linreg degree: " << degree << std::endl;

    uint64_t* const x_accum = new uint64_t[num_x + num_quad + 1];
    uint64_t* const y_accum = new uint64_t[num_x + 1];
    memset(x_accum, 0, (num_x + num_quad + 1) * sizeof(uint64_t));
    memset(y_accum, 0, (num_x + 1) * sizeof(uint64_t));

    int num_bytes = 0;
    initMsg msg;
    msg.num_of_inputs = numreqs;
    msg.type = LINREG_OP;

    if (CLIENT_BATCH) {
        num_bytes += lin_reg_helper(protocol, numreqs, degree, x_accum, y_accum, &msg);
    } else {
        auto start = clock_start();
        for (unsigned int i = 0; i < numreqs; i++)
            num_bytes += lin_reg_helper(protocol, 1, degree, x_accum, y_accum,
                                        i == 0 ? &msg : nullptr);
        std::cout << "make+send:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;
    }

    // compute answer
    double* c = SolveLinReg(degree, x_accum, y_accum);

    delete[] x_accum;
    delete[] y_accum;

    std::cout << "Estimate: y = ";
    for (unsigned int i = 0; i < degree; i++) {
        if (i > 0) std::cout << " + ";
        std::cout << c[i];
        if (i > 0) std::cout << " * x_" << (i-1);
    }
    std::cout << std::endl;

    delete[] c;

    std::cout << "Total sent bytes: " << num_bytes << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cout << "Usage: ./bin/client num_submissions server0_port server1_port OPERATION num_bits linreg_degree" << endl;
    }

    const int numreqs = atoi(argv[1]);  // Number of simulated clients
    const int port0 = atoi(argv[2]);
    const int port1 = atoi(argv[3]);

    const std::string protocol(argv[4]);

    if (argc >= 6) {
        num_bits = atoi(argv[5]);
        std::cout << "num bits: " << num_bits << std::endl;
        max_int = 1ULL << num_bits;
        std::cout << "max int: " << max_int << std::endl;
        if (num_bits > 63)
            error_exit("Num bits is too large. Int math is done mod 2^64.");
    }

    if (argc >= 7) {
        linreg_degree = atoi(argv[6]);
        std::cout << "linreg degree: " << num_bits << std::endl;
        if (linreg_degree < 2)
            error_exit("Linreg Degree must be >= 2");
    }

    // Set up server connections

    struct sockaddr_in server1, server0;

    sockfd0 = socket(AF_INET, SOCK_STREAM, 0);
    sockfd1 = socket(AF_INET, SOCK_STREAM, 0);

    if (sockfd0 < 0 or sockfd1 < 0) error_exit("Socket creation failed!");
    int sockopt = 0;
    if (setsockopt(sockfd0, SOL_SOCKET, SO_REUSEADDR, &sockopt, sizeof(sockopt)))
        error_exit("Sockopt on 0 failed");
    if (setsockopt(sockfd1, SOL_SOCKET, SO_REUSEADDR, &sockopt, sizeof(sockopt)))
        error_exit("Sockopt on 1 failed");
    if (setsockopt(sockfd0, SOL_SOCKET, SO_REUSEPORT, &sockopt, sizeof(sockopt)))
        error_exit("Sockopt on 0 failed");
    if (setsockopt(sockfd1, SOL_SOCKET, SO_REUSEPORT, &sockopt, sizeof(sockopt)))
        error_exit("Sockopt on 1 failed");

    server1.sin_port = htons(port1);
    server0.sin_port = htons(port0);

    server0.sin_family = AF_INET;
    server1.sin_family = AF_INET;

    inet_pton(AF_INET, SERVER0_IP, &server0.sin_addr);
    inet_pton(AF_INET, SERVER1_IP, &server1.sin_addr);
    std::cout << "Connecting to server 0" << std::endl;
    if (connect(sockfd0, (sockaddr*)&server0, sizeof(server0)) < 0)
        error_exit("Can't connect to server0");
    std::cout << "Connecting to server 1" << std::endl;
    if (connect(sockfd1, (sockaddr*)&server1, sizeof(server1)) < 0)
        error_exit("Can't connect to server1");

    init_constants();

    auto start = clock_start();
    if (protocol == "BITSUM") {
        std::cout << "Uploading all BITSUM shares: " << numreqs << std::endl;
        if (DEBUG_INVALID)
            bit_sum_invalid(protocol, numreqs);
        else
            bit_sum(protocol, numreqs);
        std::cout << "Total time:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;
    }

    else if (protocol == "INTSUM") {
        std::cout << "Uploading all INTSUM shares: " << numreqs << std::endl;
        if (DEBUG_INVALID)
            int_sum_invalid(protocol, numreqs);
        else
            int_sum(protocol, numreqs);
        std::cout << "Total time:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;
    }

    else if (protocol == "ANDOP") {
        std::cout << "Uploading all AND shares: " << numreqs << std::endl;
        if (DEBUG_INVALID)
            xor_op_invalid(protocol, numreqs);
        else
            xor_op(protocol, numreqs);
        std::cout << "Total time:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;
    }

    else if (protocol == "OROP") {
        std::cout << "Uploading all OR shares: " << numreqs << std::endl;
        if (DEBUG_INVALID)
            xor_op_invalid(protocol, numreqs);
        else
            xor_op(protocol, numreqs);
        std::cout << "Total time:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;
    }

    else if (protocol == "MAXOP") {
        std::cout << "Uploading all MAX shares: " << numreqs << std::endl;
        if (DEBUG_INVALID)
            max_op_invalid(protocol, numreqs);
        else
            max_op(protocol, numreqs);
        std::cout << "Total time:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;
    }

    else if (protocol == "MINOP") {
        // Min(x) = - max(-x) = b - max(b - x)
        std::cout << "Uploading all MIN shares: " << numreqs << std::endl;
        if (DEBUG_INVALID)
            max_op_invalid(protocol, numreqs);
        else
            max_op(protocol, numreqs);
        std::cout << "Total time:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;
    }

    else if (protocol == "VAROP") {
        std::cout << "Uploading all VAR shares: " << numreqs << std::endl;
        if (DEBUG_INVALID)
            var_op_invalid(protocol, numreqs);
        else
            var_op(protocol, numreqs);
        std::cout << "Total time:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;
    }

    else if (protocol == "STDDEVOP") {
        // Stddev(x) = sqrt(Var(x))
        std::cout << "Uploading all STDDEV shares: " << numreqs << std::endl;
        if (DEBUG_INVALID)
            var_op_invalid(protocol, numreqs);
        else
            var_op(protocol, numreqs);
        std::cout << "Total time:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;
    }

    else if(protocol == "LINREGOP") {
        std::cout << "Uploading all LINREG shares: " << numreqs << std::endl;

        // TODO: invalid
        lin_reg(protocol, numreqs);
        std::cout << "Total time:\t" << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;
    }

    else {
        std::cout << "Unrecognized protocol: " << protocol << std::endl;
        initMsg msg;
        msg.type = NONE_OP;
        send_to_server(0, &msg, sizeof(initMsg));
        send_to_server(1, &msg, sizeof(initMsg));
    }

    close(sockfd0);
    close(sockfd1);

    clear_constants();

    return 0;
}
