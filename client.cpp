/*
Simulates a group of num_submission clients that communicate with the servers.
*/

// TODO: Eventually htonl/ntohl wrappers on shares. Fine when client/server on same machine.

#include "client.h"

#include <arpa/inet.h>
#include <math.h>  // sqrt
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


uint64_t max_int;
uint64_t small_max_int; // sqrt(max_int)
int sockfd0, sockfd1;
bool include_invalid = false;

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

void send_maxshare(const MaxShare& maxshare, const int server_num, const unsigned int B) {
    int sock = (server_num == 0) ? sockfd0 : sockfd1;

    int s = send(sock, (void *)&maxshare, sizeof(MaxShare), 0);
    // std::cout << "sent : " << s << " bytes" << std::endl;

    for (unsigned int i = 0; i <= B; i++)
        s = send(sock, (void *)&(maxshare.arr[i]), sizeof(uint32_t), 0);
}

// Wrapper around send, with error catching.
int send_to_server(const int server, const void* const buffer, const size_t n, const int flags = 0) {
    int socket = (server == 0 ? sockfd0 : sockfd1);
    int ret = send(socket, buffer, n, flags);
    if (ret < 0) error_exit("Failed to send to server ");
    return ret;
}

void bit_sum(const std::string protocol, const size_t numreqs) {
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
    for (unsigned int i = 0; i < numreqs; i++) {
        BitShare share0, share1;
        const std::string pk_s = pub_key_to_hex((uint64_t*)&b[i]);
        const char* pk = pk_s.c_str();

        shares1[i] = real_vals[i]^shares0[i];
        ans += (real_vals[i] ? 1 : 0);

        memcpy(share0.pk, &pk[0], PK_LENGTH);
        share0.val = shares0[i];
        memcpy(share0.signature, &pk[0], PK_LENGTH);

        memcpy(share1.pk, &pk[0], PK_LENGTH);
        share1.val = shares1[i];
        memcpy(share1.signature, &pk[0], PK_LENGTH);

        send_to_server(0, &share0, sizeof(BitShare));
        send_to_server(1, &share1, sizeof(BitShare));
    }
    std::cout << "Uploaded all shares. Ans : " << ans << std::endl;

    delete[] b;
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
        const char* pk = pk_str.c_str();

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
        memcpy(share0.signature, &pk[0], PK_LENGTH);
        if (i == 0)
            share0.pk[0] = 'q';
        if (i == 2)
            memcpy(share0.pk, &prev_pk[0], PK_LENGTH);

        memcpy(share1.pk, &pk[0], PK_LENGTH);
        share1.val = shares1[i];
        memcpy(share1.signature, &pk[0], PK_LENGTH);
        if (i == 4)
            memcpy(share1.pk, &prev_pk[0], PK_LENGTH);

        send_to_server(0, &share0, sizeof(BitShare));
        send_to_server(1, &share1, sizeof(BitShare));
    }
    std::cout << "Uploaded all shares. Ans : " << ans << std::endl;

    delete[] b;
}

void int_sum(const std::string protocol, const size_t numreqs) {
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

    for (unsigned int i = 0; i < numreqs; i++) {
        real_vals[i] = real_vals[i] % max_int;
        shares0[i] = shares0[i] % max_int;
        shares1[i] = real_vals[i] ^ shares0[i];
        ans += real_vals[i];

        // std::cout << "real_vals[" << i << "] = " << real_vals[i] << " = " << shares0[i] << " ^ " << shares1[i] << std::endl;

        IntShare share0, share1;
        const std::string pk_s = pub_key_to_hex((uint64_t*)&b[i]);
        const char* pk = pk_s.c_str();

        memcpy(share0.pk, &pk[0], PK_LENGTH);
        share0.val = shares0[i];
        memcpy(share0.signature, &pk[0], PK_LENGTH);

        memcpy(share1.pk, &pk[0], PK_LENGTH);
        share1.val = shares1[i];
        memcpy(share1.signature, &pk[0], PK_LENGTH);

        send_to_server(0, &share0, sizeof(IntShare));
        send_to_server(1, &share1, sizeof(IntShare));
    }

    std::cout << "Uploaded all shares. Ans : " << ans << std::endl;

    delete[] b;
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
        const char* pk = pk_str.c_str();

        memcpy(share0.pk, &pk[0], PK_LENGTH);
        share0.val = shares0[i];
        memcpy(share0.signature, &pk[0], PK_LENGTH);
        if (i == 2)
            share0.pk[0] = 'q';
        if (i == 4)
            memcpy(share0.pk, &prev_pk[0], PK_LENGTH);

        memcpy(share1.pk, &pk[0], PK_LENGTH);
        share1.val = shares1[i];
        memcpy(share1.signature, &pk[0], PK_LENGTH);
        if (i == 6)
            memcpy(share1.pk, &prev_pk[0], PK_LENGTH);

        send_to_server(0, &share0, sizeof(IntShare));
        send_to_server(1, &share1, sizeof(IntShare));
    }

    std::cout << "Uploaded all shares. Ans : " << ans << std::endl;

    delete[] b;
}

void xor_op(const std::string protocol, const size_t numreqs) {
    initMsg msg;
    msg.num_of_inputs = numreqs;
    bool ans;
    if (protocol == "ANDOP") {
        msg.type = AND_OP;
        ans = true;
    }
    if (protocol == "OROP") {
        msg.type = OR_OP;
        ans = false;
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

    // Share splitting. Same as int sum. Sum of shares = encoded value
    for (unsigned int i = 0; i < numreqs; i++) {
        shares1[i] = encoded_values[i] ^ shares0[i];

        IntShare share0, share1;
        const std::string pk_s = pub_key_to_hex((uint64_t*)&b[i]);
        const char* pk = pk_s.c_str();

        memcpy(share0.pk, &pk[0], PK_LENGTH);
        share0.val = shares0[i];
        memcpy(share0.signature, &pk[0], PK_LENGTH);

        memcpy(share1.pk, &pk[0], PK_LENGTH);
        share1.val = shares1[i];
        memcpy(share1.signature, &pk[0], PK_LENGTH);

        send_to_server(0, &share0, sizeof(IntShare));
        send_to_server(1, &share1, sizeof(IntShare));
    }

    std::cout << "Uploaded all shares. Ans : " << std::boolalpha << ans << std::endl;

    delete[] b;
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
    }
    if (protocol == "OROP") {
        msg.type = OR_OP;
        ans = false;
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
            std::cout << " (invalid) " << std::endl;;
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
        const char* pk = pk_str.c_str();

        memcpy(share0.pk, &pk[0], PK_LENGTH);
        share0.val = shares0[i];
        memcpy(share0.signature, &pk[0], PK_LENGTH);
        if (include_invalid && i == 0)
            share0.pk[0] = 'q';
        if (include_invalid && i == 2)
            memcpy(share0.pk, &prev_pk[0], PK_LENGTH);

        memcpy(share1.pk, &pk[0], PK_LENGTH);
        share1.val = shares1[i];
        memcpy(share1.signature, &pk[0], PK_LENGTH);
        if (include_invalid && i == 4)
            memcpy(share1.pk, &prev_pk[0], PK_LENGTH);

        send_to_server(0, &share0, sizeof(IntShare));
        send_to_server(1, &share1, sizeof(IntShare));
    }

    std::cout << "Uploaded all shares. Ans : " << std::boolalpha << ans << std::endl;

    delete[] b;
}

void max_op(const std::string protocol, const size_t numreqs) {
    const unsigned int B = 250;

    initMsg msg;
    msg.num_of_inputs = numreqs;
    msg.max_inp = B;
    uint32_t ans;
    if (protocol == "MAXOP") {
        msg.type = MAX_OP;
        ans = 0;
    }
    if (protocol == "MINOP") {
        msg.type = MIN_OP;
        ans = B;
    }
    send_to_server(0, &msg, sizeof(initMsg), 0);
    send_to_server(1, &msg, sizeof(initMsg), 0);

    emp::PRG prg(fix_key);

    emp::block* const b = new block[numreqs];
    prg.random_block(b, numreqs);

    uint32_t values[numreqs];
    uint32_t or_encoded_array[B+1];
    uint32_t shares0[B+1];
    uint32_t shares1[B+1];
    prg.random_data(values, numreqs * sizeof(uint32_t));

    for (unsigned int i = 0; i < numreqs; i++) {
        MaxShare share0, share1;
        values[i] = values[i] % (B + 1);
        if (protocol == "MAXOP")
            ans = (values[i] > ans? values[i] : ans);
        if (protocol == "MINOP")
            ans = (values[i] < ans? values[i] : ans);

        prg.random_data(or_encoded_array, (B+1)*sizeof(uint32_t));
        prg.random_data(shares0, (B+1)*sizeof(uint32_t));

        // min(x) = -max(-x) = B - max(B - x)
        int v;
        if (protocol == "MAXOP")
            v = values[i];
        if (protocol == "MINOP")
            v = B - values[i];

        for (unsigned int j = v + 1; j <= B ; j++)
            or_encoded_array[j] = 0;

        for (unsigned int j = 0; j <= B; j++)
            shares1[j] = shares0[j] ^ or_encoded_array[j];

        const std::string pk_s = pub_key_to_hex((uint64_t*)&b[i]);
        const char* pk = pk_s.c_str();

        memcpy(share0.pk, &pk[0], PK_LENGTH);
        memcpy(share0.signature, &pk[0], PK_LENGTH);
        share0.arr = shares0;

        memcpy(share1.pk, &pk[0], PK_LENGTH);
        memcpy(share1.signature, &pk[0], PK_LENGTH);
        share1.arr = shares1;

        send_maxshare(share0, 0, B);
        send_maxshare(share1, 1, B);
    }

    std::cout << "Uploaded all shares. Ans : " << ans << std::endl;

    delete[] b;
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
    int ans;
    if (protocol == "MAXOP") {
        msg.type = MAX_OP;
        ans = 0;
    }
    if (protocol == "MINOP") {
        msg.type = MIN_OP;
        ans = B;
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
        int v;
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
        const char* pk = pk_str.c_str();

        memcpy(share0.pk, &pk[0], PK_LENGTH);
        memcpy(share0.signature, &pk[0], PK_LENGTH);
        share0.arr = shares0;
        if (i == 0)
            share0.pk[0] = 'q';
        if (i == 2)
            memcpy(share0.pk, &prev_pk[0], PK_LENGTH);

        memcpy(share1.pk, &pk[0], PK_LENGTH);
        memcpy(share1.signature, &pk[0], PK_LENGTH);
        share1.arr = shares1;
        if (i == 4)
            memcpy(share1.pk, &prev_pk[0], PK_LENGTH);

        send_maxshare(share0, 0, B);
        send_maxshare(share1, 1, B);
    }

    std::cout << "Uploaded all shares. Ans : " << ans << std::endl;

    delete[] b;
}

void var_op(const std::string protocol, const size_t numreqs) {
    initMsg msg;
    msg.num_of_inputs = numreqs;
    if (protocol == "VAROP")
        msg.type = VAR_OP;
    if (protocol == "STDDEVOP")
        msg.type = STDDEV_OP;
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
    uint64_t sum = 0, sumsquared = 0;

    emp::PRG prg(fix_key);
    prg.random_block(b, numreqs);
    prg.random_data(real_vals, numreqs * sizeof(uint64_t));
    prg.random_data(shares0, numreqs * sizeof(uint64_t));
    prg.random_data(shares0_squared, numreqs * sizeof(uint64_t));

    fmpz_t inp[2];
    fmpz_init(inp[0]);
    fmpz_init(inp[1]);

    for (unsigned int i = 0; i < numreqs; i++) {
        real_vals[i] = real_vals[i] % small_max_int;
        shares0[i] = shares0[i] % small_max_int;
        shares1[i] = real_vals[i] ^ shares0[i];
        uint64_t squared = real_vals[i] * real_vals[i];
        shares0_squared[i] = shares0_squared[i] % max_int;
        shares1_squared[i] = squared ^ shares0_squared[i];
        sum += real_vals[i];
        sumsquared += squared;

        // std::cout << i << ": " << real_vals[i] << " = " << shares0[i] << "^" << shares1[i];
        // std::cout << ", " << squared << " = " << shares0_squared[i] << "^" << shares1_squared[i] << std::endl;

        VarShare share0, share1;
        const std::string pk_s = pub_key_to_hex((uint64_t*)&b[i]);
        const char* pk = pk_s.c_str();

        memcpy(share0.pk, &pk[0], PK_LENGTH);
        share0.val = shares0[i];
        share0.val_squared = shares0_squared[i];
        memcpy(share0.signature, &pk[0], PK_LENGTH);

        memcpy(share1.pk, &pk[0], PK_LENGTH);
        share1.val = shares1[i];
        share1.val_squared = shares1_squared[i];
        memcpy(share1.signature, &pk[0], PK_LENGTH);

        send_to_server(0, &share0, sizeof(VarShare));
        send_to_server(1, &share1, sizeof(VarShare));
        // SNIP: proof that x^2 = x_squared
        fmpz_set_si(inp[0], real_vals[i]);
        fmpz_set_si(inp[1], real_vals[i] * real_vals[i]);
        Circuit* circuit = CheckVar();
        // Run through circuit to set wires
        circuit->Eval(inp);
        ClientPacket p0, p1;
        share_polynomials(circuit, p0, p1);
        send_ClientPacket(sockfd0, p0);
        send_ClientPacket(sockfd1, p1);
        delete p0;
        delete p1;
    }

    std::cout << "sum: " << sum << std::endl;
    std::cout << "sumsquared: " << sumsquared << std::endl;

    const double ex = 1. * sum / numreqs;
    const double ex2 = 1. * sumsquared / numreqs;
    double ans = ex2 - (ex * ex);
    std::cout << "E[X^2] - E[X]^2 = " << ex2 << " - (" << ex << ")^2 = " << ans << std::endl;
    if (protocol == "STDDEVOP")
        ans = sqrt(ans);
    std::cout << "True Ans: " << ans << std::endl;

    delete[] b;
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
    if (protocol == "VAROP")
        msg.type = VAR_OP;
    if (protocol == "STDDEVOP")
        msg.type = STDDEV_OP;
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

    std::cout << "small_max_int: " << small_max_int << std::endl;
    std::cout << "max_int: " << max_int << std::endl;

    fmpz_t inp[2];
    fmpz_init(inp[0]);
    fmpz_init(inp[1]);

    std::string pk_str = "";

    for (unsigned int i = 0; i < numreqs; i++) {
        if (i != 0)  // x not capped
            real_vals[i] = real_vals[i] % small_max_int;
        if (i != 1)  // x shares not capped
            shares0[i] = shares0[i] % small_max_int;
        shares1[i] = real_vals[i] ^ shares0[i];
        uint64_t squared = real_vals[i] * real_vals[i];
        if (i == 3 or i == 4)  // x^2 != x * x
            squared = (squared + 10) % max_int;
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
        const char* pk = pk_str.c_str();

        memcpy(share0.pk, &pk[0], PK_LENGTH);
        share0.val = shares0[i];
        share0.val_squared = shares0_squared[i];
        memcpy(share0.signature, &pk[0], PK_LENGTH);
        if (i == 9)
            share0.pk[0] = 'q';
        if (i == 11)
            memcpy(share0.pk, &prev_pk[0], PK_LENGTH);

        memcpy(share1.pk, &pk[0], PK_LENGTH);
        share1.val = shares1[i];
        share1.val_squared = shares1_squared[i];
        memcpy(share1.signature, &pk[0], PK_LENGTH);
        if (i == 13)
            memcpy(share1.pk, &prev_pk[0], PK_LENGTH);

        send_to_server(0, &share0, sizeof(VarShare));
        send_to_server(1, &share1, sizeof(VarShare));
        // SNIP: proof that x^2 = x_squared
        fmpz_set_si(inp[0], real_vals[i]);
        fmpz_set_si(inp[1], real_vals[i] * real_vals[i]);
        if (i == 3)
            fmpz_set_si(inp[1], (real_vals[i] * real_vals[i] + 10) % max_int);
        Circuit* circuit = CheckVar();
        // Run through circuit to set wires
        circuit->Eval(inp);
        ClientPacket p0, p1;
        share_polynomials(circuit, p0, p1);
        if (i == 5)
            p0->NWires = 1;  // N is used in send wrapper, so has trouble
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
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cout << "Usage: ./bin/client num_submissions server0_port server1_port OPERATION num_bits include_invalid" << endl;
    }

    const int numreqs = atoi(argv[1]);  // Number of simulated clients
    const int port0 = atoi(argv[2]);
    const int port1 = atoi(argv[3]);

    const std::string protocol(argv[4]);

    if (argc >= 6) {
        int num_bits = atoi(argv[5]);
        max_int = 1ULL << num_bits;
        small_max_int = 1ULL << (num_bits / 2);
    }

    if (argc >= 7) {
        std::stringstream ss(argv[6]);
        if (!(ss >> std::boolalpha >> include_invalid)) {
            error_exit("Could not parse to bool");
        }
        std::cout << "Include Invalid: " << std::boolalpha << include_invalid << std::endl;
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

    if (protocol == "BITSUM") {
        std::cout << "Uploading all BITSUM shares: " << numreqs << std::endl;
        if (include_invalid)
            bit_sum_invalid(protocol, numreqs);
        else
            bit_sum(protocol, numreqs);
    }

    else if (protocol == "INTSUM") {
        std::cout << "Uploading all INTSUM shares: " << numreqs << std::endl;
        if (include_invalid)
            int_sum_invalid(protocol, numreqs);
        else
            int_sum(protocol, numreqs);
    }

    else if (protocol == "ANDOP") {
        std::cout << "Uploading all AND shares: " << numreqs << std::endl;
        if (include_invalid)
            xor_op_invalid(protocol, numreqs);
        else
            xor_op(protocol, numreqs);
    }

    else if (protocol == "OROP") {
        std::cout << "Uploading all OR shares: " << numreqs << std::endl;
        if (include_invalid)
            xor_op_invalid(protocol, numreqs);
        else
            xor_op(protocol, numreqs);
    }

    else if (protocol == "MAXOP") {
        std::cout << "Uploading all MAX shares: " << numreqs << std::endl;
        if (include_invalid)
            max_op_invalid(protocol, numreqs);
        else
            max_op(protocol, numreqs);
    }

    else if (protocol == "MINOP") {
        // Min(x) = - max(-x) = b - max(b - x)
        std::cout << "Uploading all MIN shares: " << numreqs << std::endl;
        if (include_invalid)
            max_op_invalid(protocol, numreqs);
        else
            max_op(protocol, numreqs);
    }

    else if (protocol == "VAROP") {
        std::cout << "Uploading all VAR shares: " << numreqs << std::endl;
        if (include_invalid)
            var_op_invalid(protocol, numreqs);
        else
            var_op(protocol, numreqs);
    }

    else if (protocol == "STDDEVOP") {
        // Stddev(x) = sqrt(Var(x))
        std::cout << "Uploading all STDDEV shares: " << numreqs << std::endl;
        if (include_invalid)
            var_op_invalid(protocol, numreqs);
        else
            var_op(protocol, numreqs);
    }

    else if(protocol == "LINREGOP") {

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

    return 0;
}
