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


std::vector<BitShare> bitshares;
std::unordered_map<std::string,int> bitshare_map;
std::unordered_map<std::string,bool> bitshare_valid_map;
std::vector<IntShare> intshares;
std::unordered_map<std::string,uint32_t> intshare_map;
std::unordered_map<std::string,bool> intshare_valid_map;

uint32_t int_sum_max;
uint32_t num_bits;

int main(int argc, char** argv){
    if(argc < 4){
        std::cout << "Usage: ./bin/server server_num(0/1) this_port other_server_port INT_SUM_MAX_bits" << endl;
    }

    int server_num = atoi(argv[1]);
    int port = atoi(argv[2]);
    int other_port = atoi(argv[3]);

    int sockfd, newsockfd;

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

    sockaddr_in addr;

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

                bytes_read = recv(newsockfd,(void*)&bitshare,sizeof(BitShare),0);
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

                    io = new NetIO("127.0.0.1",60051);
                    uint64_t b = bitsum_ot_receiver<NetIO,OTNP>(io,&shares[0],bitshares.size());
                    std::cout << "From receiver: " << b << std::endl;
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
                bytes_read = recv(newsockfd,(void*)&pk[0],32,0);
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
            uint64_t a = bitsum_ot_sender<NetIO,OTNP>(io,&shares[0],&valid[0],num_ots);
            std::cout << "From sender: " << a<< std::endl;
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
                bytes_read = recv(newsockfd,(void*)&intshare,sizeof(IntShare),0);
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

                    io = new NetIO("127.0.0.1",60051);
                    uint64_t b = intsum_ot_receiver<NetIO,OTNP>(io,&shares[0],intshares.size(),num_bits);
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
                bytes_read = recv(newsockfd,(void*)&pk[0],32,0);
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
            uint64_t a = intsum_ot_sender<NetIO,OTNP>(io,&shares[0],&valid[0],num_ots,num_bits);
            std::cout << "From sender: " << a<< std::endl;
            long long t = time_from(start);
            std::cout << "Time taken : " << t << std::endl;
            delete[] shares;
            delete[] valid;
        }


        close(newsockfd);

    }



    return 0;
}