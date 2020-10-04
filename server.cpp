#include <sys/socket.h> 
#include <netinet/in.h> 
#include <cstdlib> 
#include <iostream> 
#include <unistd.h>

#include <vector>
#include <unordered_map>
#include <string>


#include "types.h"
#include "bitsum.h"

#define SERVER0_IP "127.0.0.1"
#define SERVER1_IP "127.0.0.1"

// #define SERVER0_IP "54.196.32.190"
// #define SERVER1_IP "54.187.97.97"

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


uint32_t int_sum_max;
uint32_t num_bits;

void bind_and_listen(sockaddr_in& addr, int& sockfd, int server_num, int port){
    sockfd = socket(AF_INET,SOCK_STREAM,0);

    if(sockfd == -1){
        std::cerr << "Socket creation failed" << std::endl;
        exit(EXIT_FAILURE);
    }
    int sockopt = 0;
    
    if(setsockopt(sockfd,SOL_SOCKET, SO_REUSEADDR,&sockopt,sizeof(sockopt))){
        std::cerr << "Sockopt failed" << std::endl;
        exit(EXIT_FAILURE);
    }

    addr.sin_family = AF_INET;
    addr.sin_addr.s_addr = INADDR_ANY;
    addr.sin_port = htons(port);

    if(bind(sockfd,(struct sockaddr*)&addr,sizeof(addr)) < 0){
        std::cerr << "Bind to port : " << port << " failed" << std::endl;
        exit(EXIT_FAILURE);
    }

    if(listen(sockfd,2) < 0){
        std::cerr << "Listen failed" << std::endl;
        exit(EXIT_FAILURE);
    }
    
}

int main(int argc, char** argv){
    if(argc < 4){
        std::cout << "Usage: ./bin/server server_num(0/1) this_port other_server_port INT_SUM_MAX_bits" << endl;
    }

    int server_num = atoi(argv[1]);
    int port = atoi(argv[2]);
    int other_port = atoi(argv[3]);

    int sockfd, newsockfd;

    sockaddr_in addr;

    bind_and_listen(addr,sockfd,server_num,port);

    while(1){
        auto addrlen = sizeof(addr);

        newsockfd = accept(sockfd,(struct sockaddr*)&addr,(socklen_t *)&addrlen);

        if(newsockfd < 0){
            std::cerr << "Connection creation failure" << std::endl;
            exit(EXIT_FAILURE);
        }
        
        initMsg msg;

        int bytes_read;

        bytes_read = recv(newsockfd,(void *)&msg,sizeof(initMsg),0);
        
        if(msg.type == BIT_SUM){
            BitShare bitshare;
            int num_ots = msg.num_of_inputs;
            bool *shares = new bool[msg.num_of_inputs];
            bool *valid = new bool[msg.num_of_inputs];
            for(int i = 0; i < msg.num_of_inputs; i++){
                bytes_read = 0;
                while(bytes_read < sizeof(BitShare))
                    bytes_read += recv(newsockfd,(char*)&bitshare+bytes_read,sizeof(BitShare)-bytes_read,0);
                std::string pk(bitshare.pk,bitshare.pk+32);
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

                    if(connect(sockfd_init,(sockaddr*)&server0,sizeof(server0)) < 0){
                        std::cerr << "Can't connect to server0" << std::endl;
                        exit(EXIT_FAILURE);
                    }

                    int t = send(sockfd_init,&msg,sizeof(initMsg),0);

                    for(int i = 0; i < bitshares.size(); i++){
                        send(sockfd_init,&bitshares[i].pk,32,0);
                        shares[i] = bitshares[i].val;
                    }
                    NetIO *io;

                    io = new NetIO(SERVER0_IP,60051);
                    uint64_t b = bitsum_ot_receiver<NetIO,SHOTExtension>(io,&shares[0],bitshares.size());
                    std::cout << "From receiver: " << b << std::endl;
                    send(sockfd_init,&b,sizeof(uint64_t),0);
                }
            }
            delete[] valid;
            delete[] shares;
        }

        if(msg.type == INIT_BIT_SUM){
            auto start = clock_start();
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
            uint64_t a = bitsum_ot_sender<NetIO,SHOTExtension>(io,&shares[0],&valid[0],num_ots);
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

        }

        if(msg.type == INT_SUM){
            IntShare intshare;
            int_sum_max = 1 << atoi(argv[4]);
            num_bits = atoi(argv[4]);
            int num_ots = msg.num_of_inputs;
            uint32_t *shares = new uint32_t[num_ots];
            bool *valid = new bool[num_ots];

            for(int i = 0; i < msg.num_of_inputs; i++){
                bytes_read = 0;
                while(bytes_read < sizeof(IntShare)){
                    bytes_read += recv(newsockfd,(char*)&intshare+bytes_read,sizeof(IntShare)-bytes_read,0);
                }
                // bytes_read = recv(newsockfd,(char*)&intshare,sizeof(IntShare),0);
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
                sleep(2);
                pid_t pid = fork();
                if(pid > 0){
                    initMsg msg;
                    msg.type = INIT_INT_SUM;
                    msg.num_of_inputs = intshares.size();
                    uint32_t *shares = new uint32_t[intshares.size()];
                    int sockfd_init = socket(AF_INET,SOCK_STREAM,0);

                    struct sockaddr_in server0;

                    server0.sin_port = htons(other_port);

                    server0.sin_family = AF_INET;

                    inet_pton(AF_INET,SERVER0_IP,&server0.sin_addr);

                    if(connect(sockfd_init,(sockaddr*)&server0,sizeof(server0)) < 0){
                        std::cerr << "Can't connect to server0" << std::endl;
                        exit(EXIT_FAILURE);
                    }

                    int t = send(sockfd_init,&msg,sizeof(initMsg),0);

                    for(int i = 0; i < intshares.size(); i++){
                        send(sockfd_init,&intshares[i].pk,32,0);
                        shares[i] = intshares[i].val;
                    }
                    NetIO *io;

                    io = new NetIO(SERVER0_IP,60051);
                    uint64_t b = intsum_ot_receiver<NetIO,SHOTExtension>(io,&shares[0],intshares.size(),num_bits);
                    send(sockfd_init,&b,sizeof(uint64_t),0);
                    std::cout << "From receiver: " << b << std::endl;
                }
            }
        }

        if(msg.type == INIT_INT_SUM){
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
            uint64_t a = intsum_ot_sender<NetIO,SHOTExtension>(io,&shares[0],&valid[0],num_ots,num_bits);
            uint64_t b;
            bytes_read = 0;
            while(bytes_read < sizeof(uint64_t))
                bytes_read += recv(newsockfd,(char*)&b+bytes_read,sizeof(uint64_t)-bytes_read,0);
            std::cout << "From sender: " << a<< std::endl;
            uint64_t aggr = a + b;
            std::cout << "Ans : " << aggr << std::endl;
            long long t = time_from(start);
            std::cout << "Time taken : " << t << std::endl;
            delete[] shares;
            delete[] valid;
        }

        if(msg.type == AND_OP){
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
                    uint32_t *shares = new uint32_t[andshares.size()];
                    int sockfd_init = socket(AF_INET,SOCK_STREAM,0);

                    struct sockaddr_in server0;

                    server0.sin_port = htons(other_port);

                    server0.sin_family = AF_INET;

                    inet_pton(AF_INET,SERVER0_IP,&server0.sin_addr);

                    if(connect(sockfd_init,(sockaddr*)&server0,sizeof(server0)) < 0){
                        std::cerr << "Can't connect to server0" << std::endl;
                        exit(EXIT_FAILURE);
                    }

                    int t = send(sockfd_init,&msg,sizeof(initMsg),0);

                    for(int i = 0; i < andshares.size(); i++){
                        send(sockfd_init,&andshares[i].pk,32,0);
                        shares[i] = andshares[i].val;
                    }
                    NetIO *io;

                    io = new NetIO(SERVER0_IP,60051);
                    uint64_t b = intsum_ot_receiver<NetIO,SHOTExtension>(io,&shares[0],andshares.size(),31);
                    send(sockfd_init,&b,sizeof(uint64_t),0);
                    std::cout << "From receiver: " << b << std::endl;
                }
            }
        }

        if(msg.type == INIT_AND_OP){
            std::cout << "Received INIT_AND_OP" << std::endl;
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

                if(andshare_map.find(pk_str) == andshare_map.end()){
                    valid[i] = false;
                }
                else{
                    valid[i] = true;
                    shares[i] = andshare_map[pk_str];
                }
            }

            NetIO *io;
            
            io = new NetIO(nullptr,60051);
            uint64_t a = intsum_ot_sender<NetIO,SHOTExtension>(io,&shares[0],&valid[0],num_ots,31);
            uint64_t b;
            bytes_read = 0;
            while(bytes_read < sizeof(uint64_t))
                bytes_read += recv(newsockfd,(char*)&b+bytes_read,sizeof(uint64_t)-bytes_read,0);
            std::cout << "From sender: " << a<< std::endl;
            uint64_t aggr = a + b;
            std::cout << "Ans of AND op : " << std::boolalpha << (aggr == 0) << std::endl;
            long long t = time_from(start);
            std::cout << "Time taken : " << t << std::endl;
            delete[] shares;
            delete[] valid;
        }
        
        if(msg.type == OR_OP){
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
                    uint32_t *shares = new uint32_t[orshares.size()];
                    int sockfd_init = socket(AF_INET,SOCK_STREAM,0);

                    struct sockaddr_in server0;

                    server0.sin_port = htons(other_port);

                    server0.sin_family = AF_INET;

                    inet_pton(AF_INET,SERVER0_IP,&server0.sin_addr);

                    if(connect(sockfd_init,(sockaddr*)&server0,sizeof(server0)) < 0){
                        std::cerr << "Can't connect to server0" << std::endl;
                        exit(EXIT_FAILURE);
                    }

                    int t = send(sockfd_init,&msg,sizeof(initMsg),0);

                    for(int i = 0; i < orshares.size(); i++){
                        send(sockfd_init,&orshares[i].pk,32,0);
                        shares[i] = orshares[i].val;
                    }
                    NetIO *io;

                    io = new NetIO(SERVER0_IP,60051);
                    uint64_t b = intsum_ot_receiver<NetIO,SHOTExtension>(io,&shares[0],orshares.size(),31);
                    send(sockfd_init,&b,sizeof(uint64_t),0);
                    std::cout << "From receiver: " << b << std::endl;
                }
            }
        }

        if(msg.type == INIT_OR_OP){
            std::cout << "Received INIT_OR_OP" << std::endl;
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

                if(orshare_map.find(pk_str) == orshare_map.end()){
                    valid[i] = false;
                }
                else{
                    valid[i] = true;
                    shares[i] = orshare_map[pk_str];
                }
            }

            NetIO *io;
            
            io = new NetIO(nullptr,60051);
            uint64_t a = intsum_ot_sender<NetIO,SHOTExtension>(io,&shares[0],&valid[0],num_ots,31);
            uint64_t b;
            bytes_read = 0;
            while(bytes_read < sizeof(uint64_t))
                bytes_read += recv(newsockfd,(char*)&b+bytes_read,sizeof(uint64_t)-bytes_read,0);
            std::cout << "From sender: " << a<< std::endl;
            uint64_t aggr = a + b;
            std::cout << "Ans of OR op : " << std::boolalpha << (aggr != 0) << std::endl;
            long long t = time_from(start);
            std::cout << "Time taken : " << t << std::endl;
            delete[] shares;
            delete[] valid;
        }

        if(msg.type == MAX_OP){
            MaxShare maxshare;

            int num_inputs = msg.num_of_inputs;
            int B = msg.max_inp;
            uint32_t* shares = new uint32_t[num_inputs*(msg.max_inp + 1)];
            bool *valid = new bool[num_inputs];
            int share_sz = (msg.max_inp + 1)*sizeof(uint32_t);
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

                    if(connect(sockfd_init,(sockaddr*)&server0,sizeof(server0)) < 0){
                        std::cerr << "Can't connect to server0" << std::endl;
                        exit(EXIT_FAILURE);
                    }

                    int t = send(sockfd_init,&msg,sizeof(initMsg),0);

                    for(int i = 0; i < maxshares.size(); i++){
                        send(sockfd_init,&maxshares[i].pk,32,0);
                    }
                    NetIO *io;

                    io = new NetIO(SERVER0_IP,60051);
                    auto b = maxop_ot_receiver<NetIO,SHOTExtension>(io,&shares[0],maxshares.size(),B);
                    for(int i = 0; i <= B; i++)
                        send(sockfd_init,&b[i],sizeof(uint64_t),0);
                    // std::cout << "From receiver: " << b << std::endl;
                }
            }

        }

        if(msg.type == INIT_MAX_OP){
            std::cout << "Received INIT_MAX_OP" << std::endl;
            auto start = clock_start();
            int num_inputs = msg.num_of_inputs;
            int B = msg.max_inp;
            uint32_t *shares = new uint32_t[num_inputs * (B+1)];
            bool *valid = new bool[msg.num_of_inputs];
            int share_sz = (B + 1)*sizeof(uint32_t);

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
                }
            }

            NetIO *io;
            // std::cout << "Hello -1" << endl;

            io = new NetIO(nullptr,60051);
            auto a = maxop_ot_sender<NetIO,SHOTExtension>(io,&shares[0],&valid[0],num_inputs,B);
            // std::cout << a.size() << " " << B+1 << endl;
            uint64_t* b = new uint64_t[B+1];
            
            // std::cout << "Hello" << endl;
            for(int j = 0; j <= B; j++){
                bytes_read = 0;
                while(bytes_read < sizeof(uint64_t))
                    bytes_read += recv(newsockfd,(char*)&(b[j])+bytes_read,sizeof(uint64_t)-bytes_read,0);
            }
            // std::cout << "Hello 2" << endl;
            for(int i = B; i >= 0; i--){
                a[i] = a[i] + b[i];
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
        }

        close(newsockfd);

    }

    return 0;
}