/*
Simulates a group of num_submission clients that communicate with the servers.
*/

#include <sys/socket.h> 
#include <netinet/in.h> 
#include <arpa/inet.h>
#include <cstdlib> 
#include <iostream> 
#include <unistd.h>
#include <string>
#include <sstream>

#include "types.h"
#include "proto.h"
#include "circuit.h"
#include "client.h"
#include "net_share.h"

// #define SERVER0_IP "52.87.230.64"
// #define SERVER1_IP "54.213.189.18"

#define SERVER0_IP "127.0.0.1"
#define SERVER1_IP "127.0.0.1"


uint64_t max_int;
uint64_t small_max_int; // sqrt(max_int)
int sockfd0, sockfd1;

void error_exit(const char *msg) {
    perror(msg);
    exit(EXIT_FAILURE);
}

std::string pub_key_to_hex(const uint64_t *key) {
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(16) << std::hex << key[0];
    ss << std::setfill('0') << std::setw(16) << std::hex << key[1];
    return ss.str();
}

void send_maxshare(const MaxShare& maxshare, const int server_num, const int B) {
    int sock = (server_num == 0) ? sockfd0 : sockfd1;

    int s = send(sock,(void *)&maxshare,sizeof(MaxShare),0);
    // std::cout << "sent : " << s << " bytes" << std::endl;

    for (int i = 0; i <= B; i++)
        s = send(sock,(void *)&(maxshare.arr[i]),sizeof(uint32_t),0);
}

// Wrapper around send, with error catching.
int send_to_server(const int server, const void* buffer, const size_t n, const int flags = 0) {
    int socket = (server == 0 ? sockfd0 : sockfd1);
    int ret = send(socket, buffer, n, flags);
    if (ret < 0) error_exit("Failed to send to server ");
    return ret;
}

void xor_op(const std::string protocol, const size_t numreqs) {
    initMsg msg;
    msg.num_of_inputs = numreqs;
    emp::block* b = new block[numreqs];
    uint32_t values[numreqs];
    uint32_t encoded_values[numreqs];
    uint32_t shares0[numreqs];
    uint32_t shares1[numreqs];
    uint32_t ans;
    if (protocol == "ANDOP") {
        msg.type = AND_OP;
        ans = 1;
    }
    if (protocol == "OROP") {
        msg.type = OR_OP;
        ans = 0;
    }

    emp::PRG prg(fix_key);

    prg.random_block(b, numreqs);
    prg.random_data(values, numreqs*sizeof(uint32_t));
    prg.random_data(encoded_values, numreqs*sizeof(uint32_t));
    prg.random_data(shares0, numreqs*sizeof(uint32_t));

    for (int i = 0; i < numreqs; i++) {
        values[i] = values[i]&1; // Take the parity as the real value
        // values[i] = 1;
        std::cout << "val[" << i << "] = " << values[i] << std::endl;
        if (protocol == "ANDOP")
            ans &= values[i];
        if (protocol == "OROP")
            ans |= values[i];
    }

    // encode step. set to all 0's for values that don't force the ans.
    if (protocol == "ANDOP")
        for (int i = 0; i < numreqs; i++)
            if (values[i] == 1)
                encoded_values[i] = 0;
    if (protocol == "OROP")
        for (int i = 0; i < numreqs; i++)
            if (values[i] == 0)
                encoded_values[i] = 0;

    // Share splitting. Same as int sum. Sum of shares = encoded value
    for (int i = 0; i < numreqs; i++)
        shares1[i] = encoded_values[i]^shares0[i];

    std::cerr << "NUM REQS " << numreqs << std::endl;

    send_to_server(0, &msg, sizeof(initMsg));
    send_to_server(1, &msg, sizeof(initMsg));

    for (int i = 0; i < numreqs; i++) {
        IntShare share0, share1;
        const char* pk = pub_key_to_hex((uint64_t*)&b[i]).c_str();

        memcpy(share0.pk, &pk[0], PK_LENGTH);
        share0.val = shares0[i];
        memcpy(share0.signature, &pk[0], PK_LENGTH);

        memcpy(share1.pk, &pk[0], PK_LENGTH);
        share1.val = shares1[i];
        memcpy(share1.signature, &pk[0], PK_LENGTH);

        send_to_server(0, &share0, sizeof(IntShare));
        send_to_server(1, &share1, sizeof(IntShare));
    }

    std::cout << "Uploaded all shares. " << "Ans : " << ans << std::endl;
    close(sockfd0);
    close(sockfd1);

    delete[] b;
}

void max_op(const std::string protocol, const size_t numreqs) {
    initMsg msg;
    msg.num_of_inputs = numreqs;
    msg.max_inp = 250;
    emp::PRG prg(fix_key);
    int B = msg.max_inp;
    int ans;
    if (protocol == "MAXOP") {
        msg.type = MAX_OP;
        ans = 0;
    }
    if (protocol == "MINOP") {
        msg.type = MIN_OP;
        ans = B;
    }

    std::cerr << "NUM REQS " << numreqs << std::endl;

    emp::block *b = new block[numreqs];
    prg.random_block(b,numreqs);

    uint32_t values[numreqs];
    uint32_t shares0[B+1];
    uint32_t shares1[B+1];
    uint32_t or_encoded_array[B+1];
    prg.random_data(values, numreqs*sizeof(uint32_t));

    for (int i = 0; i <= B; i++)
        or_encoded_array[i] = 0;

    send_to_server(0, &msg,sizeof(initMsg), 0);
    send_to_server(1, &msg,sizeof(initMsg), 0);

    for (int i = 0; i < numreqs; i++) {
        MaxShare maxshare0, maxshare1;
        values[i] = values[i] % (B + 1);
        std::cout << "value[" << i << "] = " << values[i] << std::endl;
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

        for (int j = v + 1; j <= B ; j++)
            or_encoded_array[j] = 0;

        for (int j = 0; j <= B; j++)
            shares1[j] = shares0[j] ^ or_encoded_array[j];

        const char* pk = pub_key_to_hex((uint64_t*)&b[i]).c_str();

        memcpy(maxshare0.pk, &pk[0], PK_LENGTH);
        memcpy(maxshare0.signature, &pk[0], PK_LENGTH);
        maxshare0.arr = shares0;

        memcpy(maxshare1.pk, &pk[0], PK_LENGTH);
        memcpy(maxshare1.signature, &pk[0], PK_LENGTH);
        maxshare1.arr = shares1;

        send_maxshare(maxshare0, 0, B);
        send_maxshare(maxshare1, 1, B);
    }

    std::cout << std::endl;
    std::cout << "Uploaded all shares. " << "Ans : " << ans << std::endl;
    close(sockfd0);
    close(sockfd1);

    delete[] b;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cout << "Usage: ./bin/client num_submissions server0_port server1_port OPERATION num_bits" << endl;
    }

    const int numreqs = atoi(argv[1]);  // Number of simulated clients
    const int port0 = atoi(argv[2]);
    const int port1 = atoi(argv[3]);

    std::string protocol(argv[4]);

    if (protocol == "INTSUM" or protocol == "VAROP") {
        int num_bits = atoi(argv[5]);
        max_int = 1 << num_bits;
        small_max_int = 1 << (num_bits / 2);
    }

    // Set up server connections

    struct sockaddr_in server1, server0;

    sockfd0 = socket(AF_INET,SOCK_STREAM,0);
    sockfd1 = socket(AF_INET,SOCK_STREAM,0);

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

    inet_pton(AF_INET,SERVER0_IP,&server0.sin_addr);
    inet_pton(AF_INET,SERVER1_IP,&server1.sin_addr);
    std::cout << "Connecting to server 0" << std::endl;
    if (connect(sockfd0, (sockaddr*)&server0, sizeof(server0)) < 0)
        error_exit("Can't connect to server0");
    std::cout << "Connecting to server 1" << std::endl;
    if (connect(sockfd1, (sockaddr*)&server1, sizeof(server1)) < 0)
        error_exit("Can't connect to server1");

    init_constants();

    if (protocol == "BITSUM") {
        emp::block *b = new block[numreqs];  // public keys
        bool shares0[numreqs];  // shares going to 0
        bool shares1[numreqs];  // shares going to 1
        bool real_vals[numreqs];  // actual values

        emp::PRG prg(fix_key);  // use fixed key
        prg.random_block(b,numreqs);  // randomize b per client
        prg.random_bool(shares0,numreqs);  // randomize server 0 shares per client
        prg.random_bool(shares1,numreqs);  // randomize server 1 shares per client

        // send initialize message to servers
        initMsg msg;
        msg.type = BIT_SUM;
        msg.num_of_inputs = numreqs;
        std::cerr << "NUM REQS " << numreqs << std::endl;
        send_to_server(0, &msg, sizeof(initMsg));
        send_to_server(1, &msg, sizeof(initMsg));

        int ans = 0;  // true answer
        for (int i = 0; i < numreqs; i++) {
            BitShare bitshare0,bitshare1;
            const char* pk = pub_key_to_hex((uint64_t*)&b[i]).c_str();

            // Client i sends share to server 0
            memcpy(bitshare0.pk, &pk[0],PK_LENGTH);
            bitshare0.val = shares0[i];
            memcpy(bitshare0.signature, &pk[0],PK_LENGTH);  // sign with pk?

            // Client i sends share to server 1
            memcpy(bitshare1.pk, &pk[0], PK_LENGTH);
            bitshare1.val = shares1[i];
            memcpy(bitshare1.signature, &pk[0], PK_LENGTH);

            send_to_server(0, &bitshare0, sizeof(BitShare));
            send_to_server(1, &bitshare1, sizeof(BitShare));

            // update truth
            real_vals[i] = shares0[i]^shares1[i];
            ans += (real_vals[i] ? 1 : 0);
        }

        std::cout << "Uploaded all shares. " << "Ans : " << ans << std::endl;
        close(sockfd0);
        close(sockfd1);

        delete[] b;
    }

    else if (protocol == "INTSUM") {
        emp::block *b = new block[numreqs];
        uint32_t shares0[numreqs];
        uint32_t shares1[numreqs];
        uint32_t real_vals[numreqs];
        uint64_t ans = 0;

        emp::PRG prg(fix_key);

        prg.random_block(b,numreqs);
        prg.random_data(shares0,numreqs*sizeof(uint32_t));
        prg.random_data(shares1,numreqs*sizeof(uint32_t));

        for (int i = 0; i < numreqs; i++) {
            shares0[i] = shares0[i]%max_int;
            shares1[i] = shares1[i]%max_int;
            real_vals[i] = shares0[i]^shares1[i];
            std::cout << "real_vals[" << i << "] = " << real_vals[i] << " = " << shares0[i] << " ^ " << shares1[i] << std::endl;
            ans += real_vals[i];
        }

        initMsg msg;
        msg.num_of_inputs = numreqs;
        msg.type = INT_SUM;

        send_to_server(0, &msg, sizeof(initMsg));
        send_to_server(1, &msg, sizeof(initMsg));

        for (int i = 0; i < numreqs; i++) {
            IntShare intshare0,intshare1;
            const char* pk = pub_key_to_hex((uint64_t*)&b[i]).c_str();

            memcpy(intshare0.pk, &pk[0], PK_LENGTH);
            intshare0.val = shares0[i];
            memcpy(intshare0.signature, &pk[0], PK_LENGTH);

            memcpy(intshare1.pk, &pk[0], PK_LENGTH);
            intshare1.val = shares1[i];
            memcpy(intshare1.signature, &pk[0], PK_LENGTH);
            std::cout << "key[" << i << "] = " << pub_key_to_hex((uint64_t*)&b[i]) << endl;

            send_to_server(0, &intshare0, sizeof(intshare0));
            send_to_server(1, &intshare1, sizeof(intshare1));
        }

        std::cout << "Uploaded all shares. True Ans : " << ans << std::endl;
        close(sockfd0);
        close(sockfd1);

        delete[] b;
    }

    else if (protocol == "ANDOP") {
        std::cout << "Uploading all AND shares: " << numreqs << std::endl;
        xor_op(protocol, numreqs);
    }

    else if (protocol == "OROP") {
        std::cout << "Uploading all OR shares: " << numreqs << std::endl;
        xor_op(protocol, numreqs);
    }

    else if (protocol == "MAXOP") {
        std::cout << "Uploading all MAX shares: " << numreqs << std::endl;
        max_op(protocol, numreqs);
    }

    else if (protocol == "MINOP") {
        // Min(x) = - max(-x) = b - max(b - x)
        std::cout << "Uploading all MIN shares: " << numreqs << std::endl;
        max_op(protocol, numreqs);
    }

    else if (protocol == "VAROP") {
        // Alternate idea: make each of these a function. Var calls intsum twice, then snip. 
        emp::block *b = new block[numreqs];
        // shares of x
        uint32_t shares0[numreqs];
        uint32_t shares1[numreqs];
        // shares of x^2
        uint32_t shares0_squared[numreqs];
        uint32_t shares1_squared[numreqs];
        uint32_t real_vals[numreqs];
        int sum = 0, sumsquared = 0;
        float ans, ex, ex2;

        emp::PRG prg(fix_key);

        prg.random_block(b,numreqs);
        prg.random_data(shares0,numreqs*sizeof(uint32_t));
        prg.random_data(shares1,numreqs*sizeof(uint32_t));
        prg.random_data(shares0_squared,numreqs*sizeof(uint32_t));

        // Since we are squaring, we have base values half as many bits.

        for (int i = 0; i < numreqs; i++) {
            shares0[i] = shares0[i] % small_max_int;
            shares1[i] = shares1[i] % small_max_int;
            real_vals[i] = shares0[i]^shares1[i];
            std::cout << "real_vals[" << i << "] = " << real_vals[i] << " = " << shares0[i] << " ^ " << shares1[i] << std::endl;

            shares0_squared[i] = shares0_squared[i] % max_int;
            uint32_t squared = real_vals[i] * real_vals[i];
            shares1_squared[i] = squared ^ shares0_squared[i];
            std::cout << "  real_vals[" << i << "]^2 = " << squared << " = " << shares0_squared[i] << " ^ " << shares1_squared[i] << std::endl;
            sum += real_vals[i];
            sumsquared += squared;
        }

        initMsg msg;
        msg.num_of_inputs = numreqs;
        msg.type = VAR_OP;

        send_to_server(0, &msg, sizeof(initMsg));
        send_to_server(1, &msg, sizeof(initMsg));

        std::cout << "numreqs: " << numreqs << std::endl;

        fmpz_t inp[2];
        fmpz_init(inp[0]);
        fmpz_init(inp[1]);

        for (int i = 0; i < numreqs; i++) {
            VarShare varshare0, varshare1;
            const char* pk = pub_key_to_hex((uint64_t*)&b[i]).c_str();

            memcpy(varshare0.pk, &pk[0], PK_LENGTH);
            varshare0.val = shares0[i];
            varshare0.val_squared = shares0_squared[i];
            memcpy(varshare0.signature, &pk[0], PK_LENGTH);

            memcpy(varshare1.pk, &pk[0], PK_LENGTH);
            varshare1.val = shares1[i];
            varshare1.val_squared = shares1_squared[i];
            memcpy(varshare1.signature, &pk[0], PK_LENGTH);

            // std::cout << "key[" << i << "] = " << pub_key_to_hex((uint64_t*)&b[i]) << endl;

            // std::cout << "  0: (" << shares0[i] << ", " << shares0_squared[i] << ")" << std::endl;
            // std::cout << "  1: (" << shares1[i] << ", " << shares1_squared[i] << ")" << std::endl;

            send_to_server(0, &varshare0, sizeof(varshare0));
            send_to_server(1, &varshare1, sizeof(varshare1));

            // SNIP: proof that x^2 = x_squared
            fmpz_set_si(inp[0], real_vals[i]);
            fmpz_set_si(inp[1], real_vals[i] * real_vals[i]);
            Circuit* var_circuit = CheckVar();
            // Run through circuit to set wires.
            var_circuit->Eval(inp);
            ClientPacket p0, p1;
            share_polynomials(var_circuit, p0, p1);

            // Send p0 to server0, p1 to server1
            send_ClientPacket(sockfd0, p0);
            send_ClientPacket(sockfd1, p1);
            delete p0;
            delete p1;
        }

        ex = 1. * sum / numreqs;
        ex2 = 1. * sumsquared / numreqs;
        ans = ex2 - (ex * ex);
        std::cout << "Uploaded all shares. True Ans : " << ex2 << " - (" << ex << ")^2 = " << ans << std::endl;
        close(sockfd0);
        close(sockfd1);

        delete[] b;
    }

    else if(protocol == "LINREGOP") {
        
    }

    else {
        std::cout << "Unrecognized protocol: " << protocol << std::endl;
    }

    return 0;
}
