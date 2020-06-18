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

#include "bitsum.h"

#define SERVER0_IP "52.91.132.239"
#define SERVER1_IP "34.209.244.29"

// #define SERVER0_IP "127.0.0.1"
// #define SERVER1_IP "127.0.0.1"

std::string pub_key_to_hex(uint64_t *key){
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(16) << std::hex << key[0];
    ss << std::setfill('0') << std::setw(16) << std::hex << key[1];
    return ss.str();
}

uint32_t max_int;

int main(int argc, char** argv){
    if(argc < 4){
        std::cout << "Usage: ./bin/client num_submissions server0_port server1_port BITSUM/INTSUM num_bits" << endl;
    }
    int numreqs = atoi(argv[1]);
    int port0 = atoi(argv[2]);
    int port1 = atoi(argv[3]);

    std::string protocol(argv[4]);

    if(protocol == "INTSUM")
        max_int = 1 << atoi(argv[5]);

    struct sockaddr_in server1, server0;

    int sockfd0,sockfd1 ;

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
        emp::block *b = new block[numreqs];
        bool* shares0 = new bool[numreqs];
        bool* shares1 = new bool[numreqs];
        bool* real_vals = new bool[numreqs];

        emp::PRG prg(fix_key);
        prg.random_block(b,numreqs);
        prg.random_bool(shares0,numreqs);
        prg.random_bool(shares1,numreqs);

        initMsg msg;
        msg.type = BIT_SUM;
        msg.num_of_inputs = numreqs;
        std::cerr << "NUM REQS " << numreqs << std::endl;
        int t = send(sockfd0,&msg,sizeof(initMsg),0);
        send(sockfd1,&msg,sizeof(msg),0);
        int ans = 0;
        for(int i = 0; i < numreqs; i++){
            BitShare bitshare0,bitshare1;
            std::string pk_str = pub_key_to_hex((uint64_t*)&b[i]);
            memcpy(bitshare0.pk,&pk_str.c_str()[0],32);
            bitshare0.val = shares0[i];
            memcpy(bitshare0.signature,&pk_str.c_str()[0],32);

            memcpy(bitshare1.pk,&pk_str.c_str()[0],32);
            bitshare1.val = shares1[i];
            memcpy(bitshare1.signature,&pk_str.c_str()[0],32);

            int s1 = send(sockfd0,(void *)&bitshare0,sizeof(bitshare0),0);
            int s2 = send(sockfd1,(void *)&bitshare1,sizeof(bitshare1),0);
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
            ans += real_vals[i];
        }

        initMsg msg;
        msg.num_of_inputs = numreqs;
        msg.type = INT_SUM;

        std::cerr << "NUM REQS " << numreqs << std::endl;

        send(sockfd0,&msg,sizeof(initMsg),0);
        send(sockfd1,&msg,sizeof(msg),0);

        for(int i = 0; i < numreqs; i++){
            IntShare intshare0,intshare1;
            memcpy(intshare0.pk,&pub_key_to_hex((uint64_t*)&b[i]).c_str()[0],32);
            intshare0.val = shares0[i];
            memcpy(intshare0.signature,&pub_key_to_hex((uint64_t*)&b[i]).c_str()[0],32);

            memcpy(intshare1.pk,&pub_key_to_hex((uint64_t*)&b[i]).c_str()[0],32);
            intshare1.val = shares1[i];
            memcpy(intshare1.signature,&pub_key_to_hex((uint64_t*)&b[i]).c_str()[0],32);
            std::cout << pub_key_to_hex((uint64_t*)&b[i]) << endl;

            int s1 = send(sockfd0,(void *)&intshare0,sizeof(intshare0),0);
            int s2 = send(sockfd1,(void *)&intshare1,sizeof(intshare1),0);
        }

        std::cout << "Uploaded all shares. " << "Ans : " << ans << std::endl;
        close(sockfd0);
        close(sockfd1);

        delete[] shares0;
        delete[] shares1;
        delete[] real_vals;
        delete[] b;

    }
    return 0;
}