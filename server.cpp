#include <sys/socket.h> 
#include <netinet/in.h> 
#include <cstdlib> 
#include <iostream> 
#include <unistd.h>

#include <vector>
#include <unordered_map>
#include <string>

#include "types.h"
#include "proto.h"
#include "server.h"
#include "net_share.h"

#define SERVER0_IP "127.0.0.1"
#define SERVER1_IP "127.0.0.1"

// #define SERVER0_IP "52.87.230.64"
// #define SERVER1_IP "54.213.189.18"

/*
shares: list of shares
share_map: map of pk -> val
share_valid_map: map of pk -> validity
*/

std::vector<BitShare> bitshares;
std::unordered_map<std::string,int> bitshare_map;
std::unordered_map<std::string,bool> bitshare_valid_map;

std::vector<IntShare> intshares;
std::unordered_map<std::string,uint32_t> intshare_map;
std::unordered_map<std::string,bool> intshare_valid_map;

std::vector<AndShare> andshares;
std::unordered_map<std::string,uint32_t> andshare_map;
std::unordered_map<std::string,bool> andshare_valid_map;

std::vector<OrShare> orshares;
std::unordered_map<std::string,uint32_t> orshare_map;
std::unordered_map<std::string,bool> orshare_valid_map;

std::vector<MaxShare> maxshares;
std::unordered_map<std::string,uint32_t*> maxshare_map;
std::unordered_map<std::string,bool> maxshare_valid_map;

std::vector<VarShare> varshares;
std::unordered_map<std::string, uint32_t> varshare_map;
std::unordered_map<std::string, uint32_t> varshare_map_squared;
std::unordered_map<std::string, bool> varshare_valid_map;

uint32_t int_sum_max;
uint32_t num_bits;

// TODO: const 60051 for netio?
// TODO: bytes_read loop helper functions? Since it's used a lot.

void error_exit(const char *msg){
    perror(msg);
    exit(EXIT_FAILURE);
}

void bind_and_listen(sockaddr_in& addr, int& sockfd, int server_num, int port){
    sockfd = socket(AF_INET,SOCK_STREAM,0);

    if(sockfd == -1)
        error_exit("Socket creation failed");
    int sockopt = 0;
    
    if(setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, &sockopt, sizeof(sockopt)))
        error_exit("Sockopt failed");
    if(setsockopt(sockfd, SOL_SOCKET, SO_REUSEPORT, &sockopt, sizeof(sockopt)))
        error_exit("Sockopt failed");

    bzero((char *) &addr, sizeof(addr));
    addr.sin_family = AF_INET;
    addr.sin_addr.s_addr = INADDR_ANY;
    addr.sin_port = htons(port);

    if(bind(sockfd,(struct sockaddr*)&addr,sizeof(addr)) < 0) {
        std::cerr << "Failed to bind to port: " << port << std::endl;
        error_exit("Bind to port failed");
    }

    if(listen(sockfd,2) < 0)
        error_exit("Listen failed");   
}

int main(int argc, char** argv){
    if(argc < 4){
        std::cout << "Usage: ./bin/server server_num(0/1) this_port other_server_port INT_SUM_MAX_bits" << endl;
    }

    int server_num = atoi(argv[1]);  // Server # 1 or # 2
    int port = atoi(argv[2]);        // port of this server
    int other_port = atoi(argv[3]);  // port of the other server

    if (argc >= 4)
        num_bits = atoi(argv[4]);

    int sockfd, newsockfd;

    sockaddr_in addr;

    bind_and_listen(addr,sockfd,server_num,port);

    init_constants();

    while(1){
        socklen_t addrlen = sizeof(addr);

        std::cout << "waiting for connection..." << std::endl;

        newsockfd = accept(sockfd,(struct sockaddr*)&addr,&addrlen);
        if(newsockfd < 0) error_exit("Connection creation failure");
        
        // Get an initMsg
        initMsg msg;
        int bytes_read;
        bytes_read = recv(newsockfd,(void *)&msg,sizeof(initMsg),0);
        
        if(msg.type == BIT_SUM){
            BitShare bitshare;  // Single share buffer
            // int num_ots = msg.num_of_inputs;  // OT per client
            bool *shares = new bool[msg.num_of_inputs];  // shares per client
            bool *valid = new bool[msg.num_of_inputs];  // validity per client
            for(int i = 0; i < msg.num_of_inputs; i++){
                // read in client's share
                bytes_read = 0;
                while(bytes_read < sizeof(BitShare))
                    bytes_read += recv(newsockfd,(char*)&bitshare+bytes_read,sizeof(BitShare)-bytes_read,0);
                std::string pk(bitshare.pk,bitshare.pk+32);
                // public key already seen, duplicate, so skip over.
                if(bitshare_map.find(pk) != bitshare_map.end())
                    continue;
                bitshares.push_back(bitshare);
                bitshare_map[pk] = bitshare.val;
                if(bitshare.val != 0 and bitshare.val != 1){
                    bitshare_valid_map[pk] = false;
                    valid[i] = false;
                }
                else{
                    bitshare_valid_map[pk] = true;
                    valid[i] = true;
                    shares[i] = bitshare.val;
                }
            }
            std::cerr << "Received " << msg.num_of_inputs << " shares" << std::endl;

            if(server_num == 1){
                pid_t pid = fork();
                if(pid > 0){
                    initMsg msg;
                    msg.type = INIT_BIT_SUM;
                    msg.num_of_inputs = bitshares.size();
                    bool *shares = new bool[msg.num_of_inputs];
                    int sockfd_init = socket(AF_INET,SOCK_STREAM,0);

                    struct sockaddr_in server0;

                    server0.sin_port = htons(other_port);

                    server0.sin_family = AF_INET;

                    inet_pton(AF_INET,SERVER0_IP,&server0.sin_addr);

                    if(connect(sockfd_init,(sockaddr*)&server0,sizeof(server0)) < 0)
                        error_exit("Can't connect to server0");

                    if (send(sockfd_init,&msg,sizeof(initMsg),0) < 0)
                        error_exit("Failed to send");

                    for(int i = 0; i < bitshares.size(); i++){
                        send(sockfd_init,&bitshares[i].pk,32,0);
                        shares[i] = bitshares[i].val;
                    }
                    NetIO *io;

                    io = new NetIO(SERVER0_IP,60051);
                    uint64_t b = bitsum_ot_receiver(io,&shares[0],bitshares.size());
                    std::cout << "From receiver: " << b << std::endl;
                    send(sockfd_init,&b,sizeof(uint64_t),0);
                }
            }
            delete[] valid;
            delete[] shares;
        } else if(msg.type == INIT_BIT_SUM){
            auto start = clock_start();  // for benchmark
            int num_ots = msg.num_of_inputs;
            bool *shares = new bool[msg.num_of_inputs];
            bool *valid = new bool[msg.num_of_inputs];
            for(int i = 0; i < msg.num_of_inputs; i++){
                char pk[32];
                bytes_read = 0;
                while(bytes_read < 32)
                    bytes_read += recv(newsockfd,(char*)&pk[0]+bytes_read,32-bytes_read,0);
                std::string pk_str(pk,pk+32);

                if(bitshare_map.find(pk_str) == bitshare_map.end()){
                    valid[i] = false;
                }
                else{
                    valid[i] = true;
                    shares[i] = bitshare_map[pk_str];
                }
            }
            NetIO *io;
            
            io = new NetIO(nullptr,60051);
            uint64_t a = bitsum_ot_sender(io,&shares[0],&valid[0],num_ots);
            std::cout << "From sender: " << a<< std::endl;
            uint64_t b;
            bytes_read = 0;
            while(bytes_read < sizeof(uint64_t))
                bytes_read += recv(newsockfd,(char*)&b+bytes_read,sizeof(uint64_t)-bytes_read,0);
            uint64_t aggr = a + b;
            std::cout << "Ans : " << aggr << std::endl;
            long long t = time_from(start);
            std::cout << "Time taken : " << t << std::endl;
            delete[] shares;
            delete[] valid;
        } else if(msg.type == INT_SUM){
            std::cout << "got INT_SUM" << std::endl;
            IntShare intshare;
            int_sum_max = 1 << atoi(argv[4]);
            int num_ots = msg.num_of_inputs;
            uint32_t *shares = new uint32_t[num_ots];
            bool *valid = new bool[num_ots];

            for(int i = 0; i < msg.num_of_inputs; i++){
                bytes_read = 0;
                while(bytes_read < sizeof(IntShare)){
                    bytes_read += recv(newsockfd,(char*)&intshare+bytes_read,sizeof(IntShare)-bytes_read,0);
                }
                std::string pk(intshare.pk,intshare.pk+32);

                if(intshare_map.find(pk) != intshare_map.end() or (intshare.val >= int_sum_max)){
                    continue; //Reject the input
                }
                else{
                    intshare_valid_map[pk] = true;
                    valid[i] = true;
                    shares[i] = intshare.val;
                }
                intshares.push_back(intshare);
                intshare_map[pk] = intshare.val;
            }

            std::cerr << "Received " << msg.num_of_inputs << " shares" << std::endl;

            if(server_num == 1){
                std::cout << "server 1 forking..." << std::endl;
                sleep(2);
                pid_t pid = fork();
                if(pid > 0){
                    std::cout << "server 1 fork to send init_int_sum" << std::endl;
                    initMsg msg;
                    msg.type = INIT_INT_SUM;
                    msg.num_of_inputs = intshares.size();
                    uint32_t *shares = new uint32_t[intshares.size()];
                    int sockfd_init = socket(AF_INET,SOCK_STREAM,0);

                    struct sockaddr_in server0;
                    server0.sin_port = htons(other_port);
                    server0.sin_family = AF_INET;

                    inet_pton(AF_INET,SERVER0_IP,&server0.sin_addr);

                    if(connect(sockfd_init,(sockaddr*)&server0,sizeof(server0)) < 0)
                        error_exit("Can't connect to server0");

                    if (send(sockfd_init,&msg,sizeof(initMsg),0) < 0)
                        error_exit("Failed to send");

                    for(int i = 0; i < intshares.size(); i++){
                        send(sockfd_init,&intshares[i].pk,32,0);
                        shares[i] = intshares[i].val;
                    }
                    NetIO *io;

                    io = new NetIO(SERVER0_IP,60051);
                    uint64_t b = intsum_ot_receiver(io,&shares[0],intshares.size(),num_bits);
                    send(sockfd_init,&b,sizeof(uint64_t),0);
                    std::cout << "Sending to server0: " << b << std::endl;
                }
            }
            delete[] shares;
            delete[] valid;
        } else if(msg.type == INIT_INT_SUM){
            std::cout << "Received INIT_INT_SUM" << std::endl;
            auto start = clock_start();
            int num_ots = msg.num_of_inputs;
            uint32_t *shares = new uint32_t[msg.num_of_inputs];
            bool *valid = new bool[msg.num_of_inputs];

            for(int i = 0; i < msg.num_of_inputs; i++){
                char pk[32];
                bytes_read = 0;
                while(bytes_read < 32)
                    bytes_read += recv(newsockfd,(char*)&pk[0]+bytes_read,32-bytes_read,0);
                std::string pk_str(pk,pk+32);

                if(intshare_map.find(pk_str) == intshare_map.end()){
                    valid[i] = false;
                }
                else{
                    valid[i] = true;
                    shares[i] = intshare_map[pk_str];
                }
            }
            NetIO *io;
            
            io = new NetIO(nullptr,60051);
            uint64_t a = intsum_ot_sender(io,&shares[0],&valid[0],num_ots,num_bits);
            uint64_t b;

            bytes_read = 0;
            while(bytes_read < sizeof(uint64_t))
                bytes_read += recv(newsockfd,(char*)&b+bytes_read,sizeof(uint64_t)-bytes_read,0);
            std::cout << "From sender: " << a << std::endl;
            std::cout << "Local sum: " << b << std::endl;
            uint64_t aggr = a + b;
            std::cout << "Ans : " << aggr << std::endl;
            long long t = time_from(start);
            std::cout << "Time taken : " << t << std::endl;
            delete[] shares;
            delete[] valid;
        } else if(msg.type == AND_OP){
            AndShare andshare;
            int num_ots = msg.num_of_inputs;

            uint32_t *shares = new uint32_t[num_ots];
            bool *valid = new bool[num_ots];

            for(int i = 0; i < num_ots; i++){
                bytes_read = 0;
                while(bytes_read < sizeof(AndShare)){
                    bytes_read += recv(newsockfd,(char*)&andshare+bytes_read,sizeof(AndShare)-bytes_read,0);
                }
                // bytes_read = recv(newsockfd,(char*)&intshare,sizeof(IntShare),0);
                std::string pk(andshare.pk,andshare.pk+32);

                
                if(andshare_map.find(pk) != andshare_map.end()){
                    continue; //Reject this input
                }
                
                andshare_valid_map[pk] = true;
                valid[i] = true;
                shares[i] = andshare.val;

                andshares.push_back(andshare);
                andshare_map[pk] = andshare.val;
            }

            std::cerr << "Received " << msg.num_of_inputs << " shares" << std::endl;

            if(server_num == 1){
                sleep(2);
                pid_t pid = fork();
                if(pid > 0){
                    initMsg msg;
                    msg.type = INIT_AND_OP;
                    msg.num_of_inputs = andshares.size();
                    // uint32_t *shares = new uint32_t[andshares.size()];
                    int sockfd_init = socket(AF_INET,SOCK_STREAM,0);

                    struct sockaddr_in server0;

                    server0.sin_port = htons(other_port);

                    server0.sin_family = AF_INET;

                    inet_pton(AF_INET,SERVER0_IP,&server0.sin_addr);

                    if(connect(sockfd_init,(sockaddr*)&server0,sizeof(server0)) < 0)
                        error_exit("Can't connect to server0");

                    if (send(sockfd_init,&msg,sizeof(initMsg),0) < 0)
                        error_exit("Failed to send");
                    uint32_t b = 0;
                    for(int i = 0; i < andshares.size(); i++){
                        send(sockfd_init,&andshares[i].pk,32,0);
                        b ^= andshares[i].val;
                    }
                    // NetIO *io;

                    // io = new NetIO(SERVER0_IP,60051);
                    // uint64_t b = intsum_ot_receiver<NetIO,SHOTExtension>(io,&shares[0],andshares.size(),31);
                    send(sockfd_init,&b,sizeof(uint32_t),0);
                    std::cout << "From receiver: " << b << std::endl;
                }
            }
        } else if(msg.type == INIT_AND_OP){
            std::cout << "Received INIT_AND_OP" << std::endl;
            auto start = clock_start();
            // int num_ots = msg.num_of_inputs;
            uint32_t *shares = new uint32_t[msg.num_of_inputs];
            bool *valid = new bool[msg.num_of_inputs];

            uint32_t a = 0;

            for(int i = 0; i < msg.num_of_inputs; i++){
                char pk[32];
                bytes_read = 0;
                while(bytes_read < 32)
                    bytes_read += recv(newsockfd,(char*)&pk[0]+bytes_read,32-bytes_read,0);
                std::string pk_str(pk,pk+32);

                if(andshare_map.find(pk_str) == andshare_map.end()){
                    valid[i] = false;
                }
                else{
                    valid[i] = true;
                    shares[i] = andshare_map[pk_str];
                    a ^= shares[i];
                }
            }

            // NetIO *io;
            
            // io = new NetIO(nullptr,60051);
            // uint64_t a = intsum_ot_sender<NetIO,SHOTExtension>(io,&shares[0],&valid[0],num_ots,31);
            uint32_t b;
            bytes_read = 0;
            while(bytes_read < sizeof(uint32_t))
                bytes_read += recv(newsockfd,(char*)&b+bytes_read,sizeof(uint32_t)-bytes_read,0);
            std::cout << "From sender: " << a<< std::endl;
            uint32_t aggr = a ^ b;
            std::cout << "Ans of AND op : " << std::boolalpha << (aggr == 0) << std::endl;
            long long t = time_from(start);
            std::cout << "Time taken : " << t << std::endl;
            delete[] shares;
            delete[] valid;
        } else if(msg.type == OR_OP){
            OrShare orshare;
            int num_ots = msg.num_of_inputs;

            uint32_t *shares = new uint32_t[num_ots];
            bool *valid = new bool[num_ots];

            for(int i = 0; i < num_ots; i++){
                bytes_read = 0;
                while(bytes_read < sizeof(OrShare)){
                    bytes_read += recv(newsockfd,(char*)&orshare+bytes_read,sizeof(OrShare)-bytes_read,0);
                }
                // bytes_read = recv(newsockfd,(char*)&intshare,sizeof(IntShare),0);
                std::string pk(orshare.pk,orshare.pk+32);

                
                if(orshare_map.find(pk) != orshare_map.end()){
                    continue; //Reject this input
                }
                
                orshare_valid_map[pk] = true;
                valid[i] = true;
                shares[i] = orshare.val;

                orshares.push_back(orshare);
                orshare_map[pk] = orshare.val;
            }

            std::cerr << "Received " << msg.num_of_inputs << " shares" << std::endl;

            if(server_num == 1){
                sleep(2);
                pid_t pid = fork();
                if(pid > 0){
                    initMsg msg;
                    msg.type = INIT_OR_OP;
                    msg.num_of_inputs = orshares.size();
                    // uint32_t *shares = new uint32_t[orshares.size()];
                    int sockfd_init = socket(AF_INET,SOCK_STREAM,0);

                    struct sockaddr_in server0;

                    server0.sin_port = htons(other_port);

                    server0.sin_family = AF_INET;

                    inet_pton(AF_INET,SERVER0_IP,&server0.sin_addr);

                    if(connect(sockfd_init,(sockaddr*)&server0,sizeof(server0)) < 0)
                        error_exit("Can't connect to server0");

                    if (send(sockfd_init,&msg,sizeof(initMsg),0) < 0)
                        error_exit("Failed to send");
                    uint32_t b = 0;
                    for(int i = 0; i < orshares.size(); i++){
                        send(sockfd_init,&orshares[i].pk,32,0);
                        b ^= orshares[i].val;
                    }
                    // NetIO *io;

                    // io = new NetIO(SERVER0_IP,60051);
                    // uint64_t b = intsum_ot_receiver<NetIO,SHOTExtension>(io,&shares[0],orshares.size(),31);
                    send(sockfd_init,&b,sizeof(uint32_t),0);
                    std::cout << "From receiver: " << b << std::endl;
                }
            }
        } else if(msg.type == INIT_OR_OP){
            std::cout << "Received INIT_OR_OP" << std::endl;
            auto start = clock_start();
            // int num_ots = msg.num_of_inputs;
            uint32_t *shares = new uint32_t[msg.num_of_inputs];
            bool *valid = new bool[msg.num_of_inputs];
            uint32_t a = 0;
            for(int i = 0; i < msg.num_of_inputs; i++){
                char pk[32];
                bytes_read = 0;
                while(bytes_read < 32)
                    bytes_read += recv(newsockfd,(char*)&pk[0]+bytes_read,32-bytes_read,0);
                std::string pk_str(pk,pk+32);

                if(orshare_map.find(pk_str) == orshare_map.end()){
                    valid[i] = false;
                }
                else{
                    valid[i] = true;
                    shares[i] = orshare_map[pk_str];
                    a ^= shares[i];
                }
            }

            // NetIO *io;
            
            // io = new NetIO(nullptr,60051);
            // uint64_t a = intsum_ot_sender<NetIO,SHOTExtension>(io,&shares[0],&valid[0],num_ots,31);
            uint32_t b;
            bytes_read = 0;
            while(bytes_read < sizeof(uint32_t))
                bytes_read += recv(newsockfd,(char*)&b+bytes_read,sizeof(uint32_t)-bytes_read,0);
            std::cout << "From sender: " << a<< std::endl;
            uint64_t aggr = a ^ b;
            std::cout << "Ans of OR op : " << std::boolalpha << (aggr != 0) << std::endl;
            long long t = time_from(start);
            std::cout << "Time taken : " << t << std::endl;
            delete[] shares;
            delete[] valid;
        } else if(msg.type == MAX_OP){
            MaxShare maxshare;

            int num_inputs = msg.num_of_inputs;
            int B = msg.max_inp;
            uint32_t* shares = new uint32_t[num_inputs*(msg.max_inp + 1)];
            bool *valid = new bool[num_inputs];
            // int share_sz = (msg.max_inp + 1)*sizeof(uint32_t);
            std::cout <<  "Num inputs : " << num_inputs << std::endl;

            for(int i = 0; i < num_inputs; i++){
                bytes_read = 0;
                // std::cout << "HELLLO MAX OP  " << i << "  " << num_inputs << std::endl;
                while(bytes_read < sizeof(MaxShare)){
                    bytes_read += recv(newsockfd,(char*)&maxshare+bytes_read,sizeof(MaxShare)-bytes_read,0);
                    // std::cout << "hell!" << std::endl;
                }

                std::string pk(maxshare.pk,maxshare.pk+32);
                bytes_read = 0;

                for(int j = 0; j <= B; j++){
                    bytes_read = 0;
                    while(bytes_read < sizeof(uint32_t))
                        bytes_read += recv(newsockfd,(char*)&(shares[i*(B+1)+j]) + bytes_read,sizeof(uint32_t)-bytes_read,0);
                    // std::cout << shares[i*(B+1)+j] << std::endl;
                }

                // while(bytes_read < share_sz){
                //     bytes_read += recv(newsockfd,(char*)&(shares[i * (B+1)]) + bytes_read,share_sz-bytes_read,0);
                // }

                std::cout << pk << " Share : " << i << std::endl; 
                // for(int j = 0; j <= B; j++)
                //     std::cout << i*(B+1) + j << "  " << shares[i*(B+1) + j] << "   ";
                // std::cout << std::endl;

                if(maxshare_valid_map.find(pk) != maxshare_valid_map.end())
                    continue;
                maxshare_valid_map[pk] = true;
                valid[i] = true;
                
                maxshare.arr = &shares[i*(B+1)];
                maxshare_map[pk] = maxshare.arr;
                maxshares.push_back(maxshare);
            }

            std::cerr << "Received " << msg.num_of_inputs << " shares" << std::endl;

            if(server_num == 1){
                sleep(2);
                pid_t pid = fork();
                if(pid > 0){
                    initMsg msg;
                    msg.type = INIT_MAX_OP;
                    msg.num_of_inputs = maxshares.size();
                    msg.max_inp = B;
                    // uint32_t *shares = new uint32_t[orshares.size()];
                    int sockfd_init = socket(AF_INET,SOCK_STREAM,0);

                    struct sockaddr_in server0;
                    server0.sin_port = htons(other_port);
                    server0.sin_family = AF_INET;

                    inet_pton(AF_INET,SERVER0_IP,&server0.sin_addr);

                    if(connect(sockfd_init,(sockaddr*)&server0,sizeof(server0)) < 0)
                        error_exit("Can't connect to server0");

                    if (send(sockfd_init,&msg,sizeof(initMsg),0) < 0)
                        error_exit("Failed to send");
                    uint32_t* b = new uint32_t[B+1];
                    for(int i = 0; i <= B; i++)
                        b[i] = 0;

                    for(int i = 0; i < maxshares.size(); i++){
                        send(sockfd_init,&maxshares[i].pk,32,0);
                        for(int j = 0; j <= B; j++)
                            b[j] ^= shares[i*(B+1) + j];
                    }
                    // NetIO *io;

                    // io = new NetIO(SERVER0_IP,60051);
                    // auto b = maxop_ot_receiver<NetIO,SHOTExtension>(io,&shares[0],maxshares.size(),B);
                    for(int i = 0; i <= B; i++)
                        send(sockfd_init,&b[i],sizeof(uint32_t),0);
                    // std::cout << "From receiver: " << b << std::endl;
                }
            }
        } else if(msg.type == INIT_MAX_OP){
            std::cout << "Received INIT_MAX_OP" << std::endl;
            auto start = clock_start();
            int num_inputs = msg.num_of_inputs;
            int B = msg.max_inp;
            uint32_t *shares = new uint32_t[num_inputs * (B+1)];
            bool *valid = new bool[msg.num_of_inputs];
            int share_sz = (B + 1)*sizeof(uint32_t);

            uint32_t* a = new uint32_t[B+1];

            for(int i = 0; i <= B; i++)
                a[i] = 0;

            for(int i = 0; i < msg.num_of_inputs; i++){
                char pk[32];
                bytes_read = 0;
                while(bytes_read < 32)
                    bytes_read += recv(newsockfd,(char*)&pk[0]+bytes_read,32-bytes_read,0);
                std::string pk_str(pk,pk+32);
                // std::cout << pk_str << std::endl;
                if(maxshare_valid_map.find(pk_str) == maxshare_valid_map.end()){
                    valid[i] = false;
                }
                else{
                    valid[i] = true;
                    
                    memcpy(shares + i*(B+1), maxshare_map[pk_str],share_sz);
                    for(int j = 0; j <= B ; j++)
                        a[j] ^= shares[i*(B+1)+j];
                }
            }

            // NetIO *io;
            // std::cout << "Hello -1" << endl;

            // io = new NetIO(nullptr,60051);
            // auto a = maxop_ot_sender<NetIO,SHOTExtension>(io,&shares[0],&valid[0],num_inputs,B);
            // std::cout << a.size() << " " << B+1 << endl;
            uint32_t* b = new uint32_t[B+1];
            
            for(int j = 0; j <= B; j++){
                bytes_read = 0;
                while(bytes_read < sizeof(uint32_t))
                    bytes_read += recv(newsockfd,(char*)&(b[j])+bytes_read,sizeof(uint32_t)-bytes_read,0);
            }
            
            for(int i = B; i >= 0; i--){
                std::cout << a[i] << " " << b[i] << std::endl;
                a[i] = a[i]^b[i];
                if(a[i] != 0){
                    std::cout << "Maximum : " << i << std::endl;
                    break;
                }
            }
            std:: cout << std::endl;
            // std::cout << "From sender: " << a<< std::endl;
            // uint64_t aggr = a + b;
            long long t = time_from(start);
            std::cout << "Time taken : " << t << std::endl;
            // delete[] shares;
            // delete[] valid;
        } else if(msg.type == VAR_OP) {
            // Alternate idea: make each of these a function. Var calls intsum twice, then snip. 
            std::cout << "VAR_OP" << std::endl;
            VarShare varshare;
            uint32_t small_max = 1 << (num_bits / 2);  // for values
            uint32_t var_max = 1 << num_bits;      // for values squared
            uint32_t num_inputs = msg.num_of_inputs;
            uint32_t shares[num_inputs];
            uint32_t shares_squared[num_inputs];
            bool valid[num_inputs];

            std::cout << "small_max: " << small_max << std::endl;
            std::cout << "var_max: " << var_max << std::endl;

            for (int i = 0; i < num_inputs; i++) {
                bytes_read = 0;
                while (bytes_read < sizeof(VarShare)) {
                    bytes_read += recv(newsockfd, (char*)&varshare + bytes_read, sizeof(VarShare) - bytes_read, 0);
                    if (bytes_read == 0) {
                        error_exit("Read 0 bytes. connection closed?");
                    }
                }
                std::string pk(varshare.pk, varshare.pk+32);

                // std::cout << "share[" << i << "] = (" << varshare.val << ", " << varshare.val_squared << ")" << std::endl;

                if((varshare_map.find(pk) != varshare_map.end())
                   or (varshare.val >= small_max) 
                   or (varshare.val_squared >= var_max)) {
                    continue;  // Reject the input
                }

                // TODO: snips?

                varshare_valid_map[pk] = true;
                valid[i] = true;
                shares[i] = varshare.val;
                shares_squared[i] = varshare.val_squared;
                varshares.push_back(varshare);
                varshare_map[pk] = varshare.val;
                varshare_map_squared[pk] = varshare.val_squared;
            }

            std::cout << "Recieved " << num_inputs << " shares" << std::endl;
            std::cout << varshares.size() << " are valid" << std::endl;

            if(server_num == 1){
                std::cout << "server 1 forking..." << std::endl;
                sleep(2);
                pid_t pid = fork();
                if(pid > 0){
                    std::cout << "server 1 fork to send INIT_VAR_OP" << std::endl;
                    initMsg msg;
                    msg.type = INIT_VAR_OP;
                    size_t num_shares = varshares.size();
                    msg.num_of_inputs = num_shares;
                    uint32_t shares[num_shares];
                    uint32_t shares_squared[num_shares];

                    int sockfd_init = socket(AF_INET, SOCK_STREAM, 0);
                    struct sockaddr_in server0;
                    server0.sin_port = htons(other_port);
                    server0.sin_family = AF_INET;
                    inet_pton(AF_INET, SERVER0_IP, &server0.sin_addr);

                    if(connect(sockfd_init,(sockaddr*)&server0,sizeof(server0)) < 0)
                        error_exit("Can't connect to server0");
                    if (send(sockfd_init, &msg, sizeof(initMsg),0) < 0)
                        error_exit("Failed to send");

                    for (int i = 0; i < num_shares; i++) {
                        send(sockfd_init, &varshares[i].pk, 32, 0);
                        shares[i] = varshares[i].val;
                        shares_squared[i] = varshares[i].val_squared;
                    }
                    NetIO* io = new NetIO(SERVER0_IP, 60051);
                    uint64_t b = intsum_ot_receiver(io, &shares[0], num_shares, num_bits);
                    std::cout << "sending b = " << b << std::endl;
                    send(sockfd_init, &b, sizeof(b), 0);

                    uint64_t b2 = intsum_ot_receiver(io, &shares_squared[0], num_shares, num_bits);
                    send(sockfd_init, &b2, sizeof(b2), 0);
                    std::cout << "sending b2 = " << b2 << std::endl;
                }
            }

            // TODO: snips

        } else if(msg.type == INIT_VAR_OP) { 
            std::cout << "Recieved INIT_VAR_OP" << std::endl;
            int num_inputs = msg.num_of_inputs;
            uint32_t shares[num_inputs];
            uint32_t shares_squared[num_inputs];
            bool valid[num_inputs];

            std::cout << "num_inputs: " << num_inputs << std::endl;

            for (int i = 0; i < num_inputs; i++) {
                char pk[32];
                bytes_read = 0;
                while (bytes_read < 32)
                    bytes_read += recv(newsockfd, (char*) &pk[0] + bytes_read, 32 - bytes_read, 0);
                std::string pk_str(pk, pk + 32);

                if(varshare_map.find(pk_str) == varshare_map.end()) {
                    valid[i] = false;
                } else {
                    valid[i] = true;
                    shares[i] = varshare_map[pk_str];
                    shares_squared[i] = varshare_map_squared[pk_str];
                }
            }

            NetIO* io = new NetIO(nullptr, 60051);
            uint64_t a = intsum_ot_sender(io, &shares[0], &valid[0], num_inputs, num_bits);
            std::cout << "got a: " << a << std::endl;
            uint64_t a2 = intsum_ot_sender(io, &shares_squared[0], &valid[0], num_inputs, num_bits);
            std::cout << "got a2: " << a2 << std::endl;
            uint64_t b, b2;

            bytes_read = 0;
            while (bytes_read < sizeof(uint64_t))
                bytes_read += recv(newsockfd, (char*)&b + bytes_read, sizeof(uint64_t) - bytes_read, 0);
            std::cout << "have b: " << b << std::endl;
            bytes_read = 0;
            while (bytes_read < sizeof(uint64_t))
                bytes_read += recv(newsockfd, (char*)&b2 + bytes_read, sizeof(uint64_t) - bytes_read, 0);
            std::cout << "have b2: " << b2 << std::endl;
            float ex = 1.0 * (a + b) / num_inputs;
            float ex2 = 1.0 * (a2 + b2) / num_inputs;
            float ans = ex2 - (ex * ex);
            std::cout << "Ans: " << ex2 << " - (" << ex << ")^2 = " << ans << std::endl;

        } else {
            std::cout << "Unrecognized message type: " << msg.type << std::endl;
        }

        std::cout << "end of loop" << std::endl << std::endl;
        close(newsockfd);

    }

    return 0;
}
