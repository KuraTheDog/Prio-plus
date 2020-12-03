/*
Simulates a group of num_submission clients that communicate with the servers.
*/

#include <sys/socket.h> 
#include <netinet/in.h> 
#include <arpa/inet.h>
#include <cstdlib> 
#include <iostream> 
#include <unistd.h>
#include <iomanip>
#include <vector>
#include <unordered_map>
#include <string>
#include <sstream>

#include "types.h"

#include "proto.h"

// #define SERVER0_IP "52.87.230.64"
// #define SERVER1_IP "54.213.189.18"

#define SERVER0_IP "127.0.0.1"
#define SERVER1_IP "127.0.0.1"


uint64_t max_int;
int sockfd0, sockfd1;

std::string pub_key_to_hex(uint64_t *key){
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(16) << std::hex << key[0];
    ss << std::setfill('0') << std::setw(16) << std::hex << key[1];
    return ss.str();
}

void send_maxshare(MaxShare& maxshare, int server_num, int B){
    int sock = (server_num == 0) ? sockfd0 : sockfd1;

    int s = send(sock,(void *)&maxshare,sizeof(MaxShare),0);
    // std::cout << "sent : " << s << " bytes" << std::endl;

    for(int i = 0; i <= B; i++)
        s = send(sock,(void *)&(maxshare.arr[i]),sizeof(uint32_t),0);
    // std::string pk_str(maxshare.pk,maxshare.pk+32);
    // std::cout << "  Send maxshare func : " <<  pk_str << "   " ;
    // for(int i = 0; i <= B; i++)
    //     std::cout << maxshare.arr[i] << "  ";
    // std::cout << endl;
    // std::cout << "sent : " << s << " bytes" << std::endl;
}

// Wrapper around send, with error catching.
void send_to_server(int server, void* buffer, size_t n, int flags) {
    int socket = (server == 0 ? sockfd0 : sockfd1);
    int ret = send(socket, buffer, n, flags);
    if(ret < 0) {
        std::cerr << "Failed to send to server " << server << std::endl;
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char** argv){
    if(argc < 4){
        std::cout << "Usage: ./bin/client num_submissions server0_port server1_port BITSUM/INTSUM num_bits" << endl;
    }

    int numreqs = atoi(argv[1]);  // Number of simulated clients
    int port0 = atoi(argv[2]);
    int port1 = atoi(argv[3]);

    std::string protocol(argv[4]);

    if(protocol == "INTSUM")
        max_int = 1 << (atoi(argv[5]));

    // Set up server connections

    struct sockaddr_in server1, server0;

    sockfd0 = socket(AF_INET,SOCK_STREAM,0);
    sockfd1 = socket(AF_INET,SOCK_STREAM,0);

    if(sockfd0 < 0 or sockfd1 < 0){
        std::cerr << "Socket creation failed!" << std::endl;
        exit(EXIT_FAILURE);
    }

    server1.sin_port = htons(port1);
    server0.sin_port = htons(port0);

    server0.sin_family = AF_INET;
    server1.sin_family = AF_INET;

    inet_pton(AF_INET,SERVER0_IP,&server0.sin_addr);
    inet_pton(AF_INET,SERVER1_IP,&server1.sin_addr);
    std::cout << "Connecting to server 0" << std::endl;
    if(connect(sockfd0,(sockaddr*)&server0,sizeof(server0)) < 0){
        std::cerr << "Can't connect to server0" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "Connecting to server 1" << std::endl;
    if(connect(sockfd1,(sockaddr*)&server1,sizeof(server1)) < 0){
        std::cerr << "Can't connect to server1" << std::endl;
        exit(EXIT_FAILURE);
    }

    if(protocol == "BITSUM"){
        emp::block *b = new block[numreqs];  // public keys
        bool* shares0 = new bool[numreqs];  // shares going to 0
        bool* shares1 = new bool[numreqs];  // shares going to 1
        bool* real_vals = new bool[numreqs];  // actual values

        emp::PRG prg(fix_key);  // use fixed key
        prg.random_block(b,numreqs);  // randomize b per client
        prg.random_bool(shares0,numreqs);  // randomize server 0 shares per client
        prg.random_bool(shares1,numreqs);  // randomize server 1 shares per client

        // send initialize message to servers
        initMsg msg;
        msg.type = BIT_SUM;
        msg.num_of_inputs = numreqs;
        std::cerr << "NUM REQS " << numreqs << std::endl;
        send_to_server(0,&msg,sizeof(msg),0);
        send_to_server(1,&msg,sizeof(msg),0);

        int ans = 0;  // true answer
        for(int i = 0; i < numreqs; i++){
            BitShare bitshare0,bitshare1;
            std::string pk_str = pub_key_to_hex((uint64_t*)&b[i]);

            // Client i sends share to server 0
            memcpy(bitshare0.pk,&pk_str.c_str()[0],32);
            bitshare0.val = shares0[i];
            memcpy(bitshare0.signature,&pk_str.c_str()[0],32);  // sign with pk?

            // Client i sends share to server 1
            memcpy(bitshare1.pk,&pk_str.c_str()[0],32);
            bitshare1.val = shares1[i];
            memcpy(bitshare1.signature,&pk_str.c_str()[0],32);

            send_to_server(0,(void *)&bitshare0,sizeof(bitshare0),0);
            send_to_server(1,(void *)&bitshare1,sizeof(bitshare1),0);

            // update truth
            real_vals[i] = shares0[i]^shares1[i];
            ans += (real_vals[i] ? 1 : 0);
        }

        std::cout << "Uploaded all shares. " << "Ans : " << ans << std::endl;
        close(sockfd0);
        close(sockfd1);


        delete[] shares0;
        delete[] shares1;
        delete[] b;
        delete[] real_vals;
    }

    if(protocol == "INTSUM"){
        emp::block *b = new block[numreqs];
        uint32_t *shares0 = new uint32_t[numreqs];
        uint32_t *shares1 = new uint32_t[numreqs];
        uint32_t *real_vals = new uint32_t[numreqs];
        uint64_t ans = 0;

        emp::PRG prg(fix_key);

        prg.random_block(b,numreqs);
        prg.random_data(shares0,numreqs*sizeof(uint32_t));
        prg.random_data(shares1,numreqs*sizeof(uint32_t));

        for(int i = 0; i < numreqs; i++){
            shares0[i] = shares0[i]%max_int;
            shares1[i] = shares1[i]%max_int;
            real_vals[i] = shares0[i]^shares1[i];
            std::cout << real_vals[i] << endl;
            ans += real_vals[i];
        }

        initMsg msg;
        msg.num_of_inputs = numreqs;
        msg.type = INT_SUM;

        std::cerr << "NUM REQS " << numreqs << std::endl;

        send_to_server(0,&msg,sizeof(msg),0);
        send_to_server(1,&msg,sizeof(msg),0);

        for(int i = 0; i < numreqs; i++){
            IntShare intshare0,intshare1;
            memcpy(intshare0.pk,&pub_key_to_hex((uint64_t*)&b[i]).c_str()[0],32);
            intshare0.val = shares0[i];
            memcpy(intshare0.signature,&pub_key_to_hex((uint64_t*)&b[i]).c_str()[0],32);

            memcpy(intshare1.pk,&pub_key_to_hex((uint64_t*)&b[i]).c_str()[0],32);
            intshare1.val = shares1[i];
            memcpy(intshare1.signature,&pub_key_to_hex((uint64_t*)&b[i]).c_str()[0],32);
            std::cout << pub_key_to_hex((uint64_t*)&b[i]) << endl;

            send_to_server(0,(void *)&intshare0,sizeof(intshare0),0);
            send_to_server(1,(void *)&intshare1,sizeof(intshare1),0);
        }

        std::cout << "Uploaded all shares. " << "Ans : " << ans << std::endl;
        close(sockfd0);
        close(sockfd1);

        delete[] shares0;
        delete[] shares1;
        delete[] real_vals;
        delete[] b;

    }

    if(protocol == "ANDOP"){
        std::cout << "Uploading all and shares. " << numreqs << std::endl;
        emp::block *b = new block[numreqs];

        uint32_t *values = new uint32_t[numreqs];
        uint32_t *encoded_values = new uint32_t[numreqs];
        uint32_t *shares0 = new uint32_t[numreqs];
        uint32_t *shares1 = new uint32_t[numreqs];
        uint64_t ans = 0;

        emp::PRG prg(fix_key);

        prg.random_block(b,numreqs);
        prg.random_data(values, numreqs*sizeof(uint32_t));
        prg.random_data(encoded_values, numreqs*sizeof(uint32_t));
        prg.random_data(shares0, numreqs*sizeof(uint32_t));

        for(int i = 0; i < numreqs; i++){
            values[i] = 1; // Take the parity as the real value
            
        }

        for(int i = 0; i < numreqs; i++){   // encode step. x_i = 1 -> enc(x_i) = 0 else enc(x_i) = r
            if(values[i] == 1){
                encoded_values[i] = 0;
            }
        }

        for(int i = 0; i < numreqs; i++)            // Share splitting. Same as int sum. Sum of shares = encoded value
            shares1[i] = encoded_values[i]^shares0[i];
        
        initMsg msg;
        msg.num_of_inputs = numreqs;
        msg.type = AND_OP;

        std::cerr << "NUM REQS " << numreqs << std::endl;

        send_to_server(0,&msg,sizeof(msg),0);
        send_to_server(1,&msg,sizeof(msg),0);

        for(int i = 0; i < numreqs; i++){
            AndShare andshare0, andshare1;

            memcpy(andshare0.pk,&pub_key_to_hex((uint64_t*)&b[i]).c_str()[0],32);
            andshare0.val = shares0[i];
            memcpy(andshare0.signature,&pub_key_to_hex((uint64_t*)&b[i]).c_str()[0],32);

            memcpy(andshare1.pk,&pub_key_to_hex((uint64_t*)&b[i]).c_str()[0],32);
            andshare1.val = shares1[i];
            memcpy(andshare1.signature,&pub_key_to_hex((uint64_t*)&b[i]).c_str()[0],32);

            send_to_server(0,(void *)&andshare0,sizeof(andshare0),0);
            send_to_server(1,(void *)&andshare1,sizeof(andshare1),0);
        }

        std::cout << "Uploaded all shares. " << "Ans : " << ans << std::endl;
        close(sockfd0);
        close(sockfd1);

        delete[] shares0;
        delete[] shares1;
        delete[] values;
        delete[] encoded_values;
        delete[] b;
    }

    if(protocol == "OROP"){
        emp::block *b = new block[numreqs];

        uint32_t *values = new uint32_t[numreqs];
        uint32_t *encoded_values = new uint32_t[numreqs];
        uint32_t *shares0 = new uint32_t[numreqs];
        uint32_t *shares1 = new uint32_t[numreqs];
        uint64_t ans = 0;

        emp::PRG prg(fix_key);

        prg.random_block(b,numreqs);
        prg.random_data(values, numreqs*sizeof(uint32_t));
        prg.random_data(encoded_values, numreqs*sizeof(uint32_t));
        prg.random_data(shares0, numreqs*sizeof(uint32_t));

        for(int i = 0; i < numreqs; i++)
            values[i] = values[i]&1; // Take the parity as the real value

        for(int i = 0; i < numreqs; i++){   // encode step. x_i = 1 -> enc(x_i) = 0 else enc(x_i) = r
            if(values[i] == 0){
                encoded_values[i] = 1;
            }
        }

        for(int i = 0; i < numreqs; i++)            // Share splitting. Same as int sum. Sum of shares = encoded value
            shares1[i] = encoded_values[i]^shares0[i];
        
        initMsg msg;
        msg.num_of_inputs = numreqs;
        msg.type = OR_OP;

        std::cerr << "NUM REQS " << numreqs << std::endl;

        send_to_server(0,&msg,sizeof(msg),0);
        send_to_server(1,&msg,sizeof(msg),0);

        for(int i = 0; i < numreqs; i++){
            OrShare orshare0, orshare1;

            memcpy(orshare0.pk,&pub_key_to_hex((uint64_t*)&b[i]).c_str()[0],32);
            orshare0.val = shares0[i];
            memcpy(orshare0.signature,&pub_key_to_hex((uint64_t*)&b[i]).c_str()[0],32);

            memcpy(orshare1.pk,&pub_key_to_hex((uint64_t*)&b[i]).c_str()[0],32);
            orshare1.val = shares1[i];
            memcpy(orshare1.signature,&pub_key_to_hex((uint64_t*)&b[i]).c_str()[0],32);

            send_to_server(0,(void *)&orshare0,sizeof(orshare0),0);
            send_to_server(1,(void *)&orshare1,sizeof(orshare1),0);
        }

        std::cout << "Uploaded all shares. " << "Ans : " << ans << std::endl;
        close(sockfd0);
        close(sockfd1);

        delete[] shares0;
        delete[] shares1;
        delete[] values;
        delete[] encoded_values;
        delete[] b;
    }

    if(protocol == "MAXOP"){
        initMsg msg;
        msg.num_of_inputs = numreqs;
        msg.type = MAX_OP;
        msg.max_inp = 250;
        emp::PRG prg(fix_key);
        int B = msg.max_inp;

        std::cerr << "NUM REQS " << numreqs << std::endl;

        emp::block *b = new block[numreqs];
        prg.random_block(b,numreqs);

        uint32_t *values = new uint32_t[numreqs];
        uint32_t *shares0 = new uint32_t[B+1];
        uint32_t *shares1 = new uint32_t[B+1];
        uint32_t *or_encoded_array = new uint32_t[B+1];
        prg.random_data(values, numreqs*sizeof(uint32_t));



        for(int i = 0; i <= msg.max_inp; i++)
            or_encoded_array[i] = 0;

        send_to_server(0,&msg,sizeof(msg),0);
        send_to_server(1,&msg,sizeof(msg),0);

        for(int i = 0; i < numreqs; i++){
            MaxShare maxshare0, maxshare1;
            values[i] = 10%(B+1);
            std::cout << values[i] << std::endl;
            prg.random_data(or_encoded_array, (B+1)*sizeof(uint32_t));
            prg.random_data(shares0, (B+1)*sizeof(uint32_t));
            for(int j = values[i] + 1; j <= B ; j++)
                or_encoded_array[j] = 0;

            for(int j = 0; j <= msg.max_inp ; j++){
                shares1[j] = shares0[j] ^ or_encoded_array[j];
                std::cout << shares0[j] << "   " << shares1[j] << endl;
            }

            memcpy(maxshare0.pk,&pub_key_to_hex((uint64_t*)&b[i]).c_str()[0],32);
            memcpy(maxshare0.signature,&pub_key_to_hex((uint64_t*)&b[i]).c_str()[0],32);
            maxshare0.arr = shares0;

            memcpy(maxshare1.pk,&pub_key_to_hex((uint64_t*)&b[i]).c_str()[0],32);
            memcpy(maxshare1.signature,&pub_key_to_hex((uint64_t*)&b[i]).c_str()[0],32);
            maxshare1.arr = shares1;

            send_maxshare(maxshare0,0,B);
            send_maxshare(maxshare1,1,B);

        }

        std::cout << std::endl;
        // std::cout << "Uploaded all shares. " << "Ans : " << ans << std::endl;
        close(sockfd0);
        close(sockfd1);

        delete[] shares0;
        delete[] shares1;
        delete[] values;
        delete[] or_encoded_array;
        delete[] b;

    }

    return 0;
}

