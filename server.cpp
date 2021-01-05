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

/* TODO: Split client and server ports
server client_in, server_in, other_server
Servers chat on server_in/other_server
So server can connect to client + other server
*/

uint32_t int_sum_max;
uint32_t num_bits;

// TODO: const 60051 for netio?

void error_exit(const char *msg) {
    perror(msg);
    exit(EXIT_FAILURE);
}

// Use example: type obj; read_in(sockfd, &obj, sizeof(type))
size_t read_in(const int sockfd, void* buf, const size_t len) {
    size_t bytes_read = 0;
    char* bufptr = (char*) buf;
    while (bytes_read < len)
        bytes_read += recv(sockfd, bufptr + bytes_read, len - bytes_read, 0);
    return bytes_read;
}

size_t send_out(const int sockfd, const void* buf, const size_t len) {
    size_t ret = send(sockfd, buf, len, 0);
    if (ret <= 0) error_exit("Failed to send");
    return ret;
}

void bind_and_listen(sockaddr_in& addr, int& sockfd, const int port, const int reuse = 1) {
    sockfd = socket(AF_INET,SOCK_STREAM,0);

    if (sockfd == -1)
        error_exit("Socket creation failed");
    
    if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, &reuse, sizeof(reuse)))
        error_exit("Sockopt failed");
    if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEPORT, &reuse, sizeof(reuse)))
        error_exit("Sockopt failed");

    bzero((char *) &addr, sizeof(addr));
    addr.sin_family = AF_INET;
    addr.sin_addr.s_addr = INADDR_ANY;
    addr.sin_port = htons(port);

    if (bind(sockfd,(struct sockaddr*)&addr,sizeof(addr)) < 0) {
        std::cerr << "Failed to bind to port: " << port << std::endl;
        error_exit("Bind to port failed");
    }

    if (listen(sockfd,2) < 0)
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
    if (sockfd == -1) error_exit("Socket creation failed");

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

bool run_snip(Circuit* circuit, const ClientPacket packet, const fmpz_t randomX, 
              const ShareSender sender, const ShareReceiver receiver, const int server_num) {
    Checker* checker = new Checker(circuit, server_num);
    checker->setReq(packet);
    CheckerPreComp* pre = new CheckerPreComp(circuit);
    pre->setCheckerPrecomp(randomX);
    CorShare* cor_share = checker->CorShareFn(pre);
    CorShare* cor_share_other = new CorShare();

    if (server_num == 0) {
        sender.CorShare(cor_share);
        receiver.CorShare(cor_share_other);
    } else {
        receiver.CorShare(cor_share_other);
        sender.CorShare(cor_share);
    }

    Cor* cor = checker->CorFn(cor_share, cor_share_other);
    fmpz_t valid_share, valid_share_other;
    fmpz_init(valid_share);
    fmpz_init(valid_share_other);
    checker->OutShare(valid_share, cor);

    if (server_num == 0) {
        sender.fmpz(valid_share);
        receiver.fmpz(valid_share_other);
    } else {
        receiver.fmpz(valid_share_other);
        sender.fmpz(valid_share);
    }

    bool circuit_valid = checker->OutputIsValid(valid_share, valid_share_other);

    fmpz_clear(valid_share);
    fmpz_clear(valid_share_other);

    return circuit_valid;
}

int main(int argc, char** argv) {
    if (argc < 5) {
        std::cout << "Usage: ./bin/server server_num(0/1) this_client_port this_server_port other_server_port INT_SUM_MAX_bits" << endl;
    }

    const int server_num = atoi(argv[1]);  // Server # 1 or # 2
    const int client_port = atoi(argv[2]); // port of this server, for the client
    const int server_port = atoi(argv[3]); // port of this server, for the other server
    const int other_port = atoi(argv[4]);  // port of the other server

    std::cout << "This server is server # " << server_num << std::endl;
    std::cout << "  Listening for client on " << client_port << std::endl;
    std::cout << "  Listening for server on " << server_port << std::endl;
    std::cout << "  Other server's port is  " << other_port  << std::endl;

    if (argc >= 5)
        num_bits = atoi(argv[5]);

    init_constants();

    // Server 0 listens, via newsockfd_server.
    // Server 1 connects, via sockfd_server
    int sockfd_server, newsockfd_server, serverfd = 0;
    if (server_num == 0) {
        server0_listen(sockfd_server, newsockfd_server, server_port, 1);
        serverfd = newsockfd_server;
    } else if (server_num == 1) {
        server1_connect(sockfd_server, other_port, 1);
        serverfd = sockfd_server;
    } else {
        error_exit("Can only handle servers #0 and #1");
    }
    ShareSender server_share_sender(serverfd);
    ShareReceiver server_share_receiver(serverfd);

    // Share the same randomX.
    // TODO: change every once in a while.
    fmpz_t randomX;
    fmpz_init(randomX);
    if (server_num == 0) {
        server_share_receiver.fmpz(randomX);
        std::cout << "Got randomX: "; fmpz_print(randomX); std::cout << std::endl;
    } else {
        fmpz_randm(randomX, seed, Int_Modulus);
        std::cout << "Sending randomX: "; fmpz_print(randomX); std::cout << std::endl;
        server_share_sender.fmpz(randomX);
    }

    int sockfd, newsockfd;
    sockaddr_in addr;

    bind_and_listen(addr, sockfd, client_port);

    while(1) {
        socklen_t addrlen = sizeof(addr);

        std::cout << "waiting for connection..." << std::endl;

        newsockfd = accept(sockfd,(struct sockaddr*)&addr,&addrlen);
        if (newsockfd < 0) error_exit("Connection creation failure");
        
        // Get an initMsg
        initMsg msg;
        read_in(newsockfd, &msg, sizeof(initMsg));
        
        if (msg.type == BIT_SUM) {
            std::cout << "BIT_SUM" << std::endl;
            auto start = clock_start();  // for benchmark

            std::vector<BitShare> bitshares;
            std::unordered_map<std::string, int> bitshare_map;

            BitShare bitshare;
            const size_t num_inputs = msg.num_of_inputs;

            for (int i = 0; i < num_inputs; i++) {
                // read in client's share
                read_in(newsockfd, &bitshare, sizeof(BitShare));
                std::string pk(bitshare.pk, bitshare.pk + 32);
                // public key already seen, duplicate, so skip over.
                if (bitshare_map.find(pk) != bitshare_map.end()
                    or (bitshare.val != 0 and bitshare.val != 1))
                    continue;
                bitshares.push_back(bitshare);
                bitshare_map[pk] = bitshare.val;
            }
            std::cerr << "Received " << msg.num_of_inputs << " shares" << std::endl;

            std::cout << "Forking for server chats" << std::endl;
            pid_t pid = fork();
            if (pid == 0) {
                std::cout << " Forked Child leaving." << std::endl;
                continue;
            }
            std::cout << " Parent for server chat" << std::endl;

            if (server_num == 1) {
                const size_t num_inputs = bitshares.size();
                send_out(serverfd, &num_inputs, sizeof(size_t));

                bool shares[num_inputs];
                for (int i = 0; i < num_inputs; i++) {
                    send_out(serverfd, &bitshares[i].pk, 32);
                    shares[i] = bitshares[i].val;
                }
                NetIO* io = new NetIO(SERVER0_IP, 60051);
                uint64_t b = bitsum_ot_receiver(io, &shares[0], num_inputs);
                std::cout << "Sending b: " << b << std::endl;
                send_out(serverfd, &b, sizeof(uint64_t));
            } else {
                size_t num_inputs;
                read_in(serverfd, &num_inputs, sizeof(size_t));
                bool shares[num_inputs];
                bool valid[num_inputs];

                for (int i = 0; i < num_inputs; i++) {
                    char pk_buf[32];
                    read_in(serverfd, &pk_buf[0], 32);
                    std::string pk(pk_buf, pk_buf + 32);

                    bool is_valid = (bitshare_map.find(pk) != bitshare_map.end());
                    valid[i] = is_valid;
                    if (!is_valid)
                        continue;
                    shares[i] = bitshare_map[pk];
                }
                NetIO* io = new NetIO(nullptr, 60051);
                uint64_t a = bitsum_ot_sender(io, &shares[0], &valid[0], num_inputs);
                std::cout << "Have a: " << a << std::endl;
                uint64_t b;
                read_in(serverfd, &b, sizeof(uint64_t));
                std::cout << "Got b: " << b << std::endl;
                uint64_t aggr = a + b;
                std::cout << "Ans : " << aggr << std::endl;
            }
            long long t = time_from(start);
            std::cout << "Time taken : " << t << std::endl;
        } else if (msg.type == INT_SUM) {
            std::cout << "INT_SUM" << std::endl;
            auto start = clock_start();  // for benchmark

            std::vector<IntShare> intshares;
            std::unordered_map<std::string, uint32_t> intshare_map;

            IntShare intshare;
            int_sum_max = 1 << num_bits;
            const size_t num_inputs = msg.num_of_inputs;

            for (int i = 0; i < num_inputs; i++) {
                read_in(newsockfd, &intshare, sizeof(IntShare));
                std::string pk(intshare.pk, intshare.pk + 32);

                if (intshare_map.find(pk) != intshare_map.end() 
                    or (intshare.val >= int_sum_max)) {
                    continue; //Reject the input
                }
                intshares.push_back(intshare);
                intshare_map[pk] = intshare.val;
            }

            std::cerr << "Received " << msg.num_of_inputs << " shares" << std::endl;

            std::cout << "Forking for server chats" << std::endl;
            pid_t pid = fork();
            if (pid == 0) {
                std::cout << " Forked Child leaving." << std::endl;
                continue;
            }
            std::cout << " Parent for server chat" << std::endl;

            if (server_num == 1) {
                const size_t num_inputs = intshares.size();
                send_out(serverfd, &num_inputs, sizeof(size_t));

                uint32_t shares[num_inputs];
                for (int i = 0; i < num_inputs; i++) {
                    send_out(serverfd, &intshares[i].pk, 32);
                    shares[i] = intshares[i].val;
                }
                NetIO* io = new NetIO(SERVER0_IP, 60051);
                uint64_t b = intsum_ot_receiver(io, &shares[0], num_inputs, num_bits);
                std::cout << "Sending b: " << b << std::endl;
                send_out(serverfd, &b, sizeof(uint64_t));
            } else {
                size_t num_inputs;
                read_in(serverfd, &num_inputs, sizeof(size_t));
                uint32_t shares[num_inputs];
                bool valid[num_inputs];

                for (int i = 0; i < num_inputs; i++) {
                    char pk_buf[32];
                    read_in(serverfd, &pk_buf[0], 32);
                    std::string pk(pk_buf, pk_buf + 32);

                    bool is_valid = (intshare_map.find(pk) != intshare_map.end());
                    valid[i] = is_valid;
                    if (!is_valid)
                        continue;
                    shares[i] = intshare_map[pk];
                }
                NetIO* io = new NetIO(nullptr, 60051);
                uint64_t a = intsum_ot_sender(io, &shares[0], &valid[0], num_inputs, num_bits);
                std::cout << "Have a: " << a << std::endl;
                uint64_t b;
                read_in(serverfd, &b, sizeof(uint64_t));
                std::cout << "Got b: " << b << std::endl;
                uint64_t aggr = a + b;
                std::cout << "Ans : " << aggr << std::endl;
            }
            long long t = time_from(start);
            std::cout << "Time taken : " << t << std::endl;
        } else if (msg.type == AND_OP) {
            std::cout << "AND_OP" << std::endl;
            auto start = clock_start();  // for benchmark

            std::vector<AndShare> andshares;
            std::unordered_map<std::string,uint32_t> andshare_map;

            AndShare andshare;
            const size_t num_inputs = msg.num_of_inputs;

            for (int i = 0; i < num_inputs; i++) {
                read_in(newsockfd, &andshare, sizeof(AndShare));
                std::string pk(andshare.pk, andshare.pk + 32);
                
                if (andshare_map.find(pk) != andshare_map.end())
                    continue;
                andshares.push_back(andshare);
                andshare_map[pk] = andshare.val;
            }

            std::cerr << "Received " << msg.num_of_inputs << " shares" << std::endl;

            std::cout << "Forking for server chats" << std::endl;
            pid_t pid = fork();
            if (pid == 0) {
                std::cout << " Forked Child leaving." << std::endl;
                continue;
            }
            std::cout << " Parent for server chat" << std::endl;

            if (server_num == 1) {
                const size_t num_inputs = andshares.size();
                send_out(serverfd, &num_inputs, sizeof(size_t));
                uint32_t b = 0;
                for (int i = 0; i < num_inputs; i++) {
                    send_out(serverfd, &andshares[i].pk, 32);
                    bool other_valid;
                    read_in(serverfd, &other_valid, sizeof(bool));
                    if (!other_valid)
                        continue;
                    b ^= andshares[i].val;
                }
                send_out(serverfd, &b, sizeof(uint32_t));
                std::cout << "Sending b: " << b << std::endl;
            } else {
                size_t num_inputs;
                read_in(serverfd, &num_inputs, sizeof(size_t));

                uint32_t a = 0;

                for (int i = 0; i < num_inputs; i++) {
                    char pk_buf[32];
                    read_in(serverfd, &pk_buf[0], 32);
                    std::string pk(pk_buf, pk_buf + 32);

                    bool is_valid = (andshare_map.find(pk) != andshare_map.end());
                    send_out(serverfd, &is_valid, sizeof(bool));
                    if (!is_valid)
                        continue;
                    a ^= andshare_map[pk];
                }
                std::cout << "Have a: " << a << std::endl;
                uint32_t b;
                read_in(serverfd, &b, sizeof(uint32_t));
                std::cout << "Got b: " << b << std::endl;
                uint32_t aggr = a ^ b;
                std::cout << "Aggr = " << aggr << std::endl;
                std::cout << "Ans: " << std::boolalpha << (aggr == 0) << std::endl;
            }
            long long t = time_from(start);
            std::cout << "Time taken : " << t << std::endl;
        } else if (msg.type == OR_OP) {
            std::cout << "OR_OP" << std::endl;
            auto start = clock_start();  // for benchmark
            std::vector<OrShare> orshares;
            std::unordered_map<std::string,uint32_t> orshare_map;

            OrShare orshare;
            const size_t num_inputs = msg.num_of_inputs;

            for (int i = 0; i < num_inputs; i++) {
                read_in(newsockfd, &orshare, sizeof(OrShare));
                std::string pk(orshare.pk, orshare.pk + 32);
                
                if (orshare_map.find(pk) != orshare_map.end())
                    continue;
                orshares.push_back(orshare);
                orshare_map[pk] = orshare.val;
            }

            std::cerr << "Received " << msg.num_of_inputs << " shares" << std::endl;

            std::cout << "Forking for server chats" << std::endl;
            pid_t pid = fork();
            if (pid == 0) {
                std::cout << " Forked Child leaving." << std::endl;
                continue;
            }
            std::cout << " Parent for server chat" << std::endl;

            if (server_num == 1) {
                const size_t num_inputs = orshares.size();
                send_out(serverfd, &num_inputs, sizeof(size_t));
                uint32_t b = 0;
                for (int i =0; i < num_inputs; i++) {
                    send_out(serverfd, &orshares[i].pk, 32);
                    bool other_valid;
                    read_in(serverfd, &other_valid, sizeof(bool));
                    if (!other_valid)
                        continue;
                    b ^= orshares[i].val;
                }
                send_out(serverfd, &b, sizeof(uint32_t));
                std::cout << "Sending b: " << b << std::endl;
            } else {
                size_t num_inputs;
                read_in(serverfd, &num_inputs, sizeof(size_t));

                uint32_t a = 0;

                for (int i = 0; i < num_inputs; i++) {
                    char pk_buf[32];
                    read_in(serverfd, &pk_buf[0], 32);
                    std::string pk(pk_buf, pk_buf + 32);

                    bool is_valid = (orshare_map.find(pk) != orshare_map.end());
                    send_out(serverfd, &is_valid, sizeof(bool));
                    if (!is_valid)
                        continue;
                    a ^= orshare_map[pk];
                }
                std::cout << "Have a: " << a << std::endl;
                uint32_t b;
                read_in(serverfd, &b, sizeof(uint32_t));
                std::cout << "Got b: " << b << std::endl;
                uint32_t aggr = a ^ b;
                std::cout << "Aggr = " << aggr << std::endl;
                std::cout << "Ans: " << std::boolalpha << (aggr != 0) << std::endl;
            }
            long long t = time_from(start);
            std::cout << "Time taken : " << t << std::endl;
        } else if (msg.type == MAX_OP) {
            std::cout << "MAX_OP" << std::endl;
            auto start = clock_start();  // for benchmark

            std::vector<MaxShare> maxshares;
            std::unordered_map<std::string, uint32_t*> maxshare_map;

            MaxShare maxshare;
            const size_t num_inputs = msg.num_of_inputs;
            const int B = msg.max_inp;
            uint32_t shares[num_inputs * (B + 1)];
            std::cout <<  "Num inputs : " << num_inputs << std::endl;

            for (int i = 0; i < num_inputs; i++) {
                read_in(newsockfd, &maxshare, sizeof(MaxShare));
                std::string pk(maxshare.pk, maxshare.pk + 32);

                for (int j = 0; j <= B; j++)
                    read_in(newsockfd, &(shares[i*(B+1)+j]), sizeof(uint32_t));


                std::cout << pk << " Share : " << i << std::endl; 
                // for(int j = 0; j <= B; j++)
                //     std::cout << i*(B+1) + j << "  " << shares[i*(B+1) + j] << "   ";
                // std::cout << std::endl;
                if (maxshare_map.find(pk) != maxshare_map.end())
                    continue;

                maxshare.arr = &shares[i*(B+1)];
                maxshare_map[pk] = maxshare.arr;
                maxshares.push_back(maxshare);
            }

            std::cerr << "Received " << msg.num_of_inputs << " shares" << std::endl;

            std::cout << "Forking for server chats" << std::endl;
            pid_t pid = fork();
            if (pid == 0) {
                std::cout << " Forked Child leaving." << std::endl;
                continue;
            }
            std::cout << " Parent for server chat" << std::endl;

            if (server_num == 1) {
                const size_t num_inputs = maxshares.size();
                send_out(serverfd, &num_inputs, sizeof(size_t));
                uint32_t b[B+1];
                for (int j = 0; j <= B; j++)
                    b[j] = 0;

                for (int i = 0; i < num_inputs; i++) {
                    send_out(serverfd, &maxshares[i].pk, 32);
                    bool other_valid;
                    read_in(serverfd, &other_valid, sizeof(bool));
                    if (!other_valid)
                        continue;
                    for (int j = 0; j <= B; j++)
                        b[j] ^= shares[i*(B+1)+j];
                }
                for (int j = 0; j <= B; j++)
                    send_out(serverfd, &b[j], sizeof(uint32_t));
            } else {
                size_t num_inputs;
                read_in(serverfd, &num_inputs, sizeof(size_t));
                uint32_t share[B+1];
                int share_sz = (B+1)*sizeof(uint32_t);
                uint32_t a[B+1];

                for (int j = 0; j <= B; j++)
                    a[j] = 0;

                for (int i = 0; i < num_inputs; i++) {
                    char pk_buf[32];
                    read_in(serverfd, &pk_buf[0], 32);
                    std::string pk(pk_buf, pk_buf + 32);

                    bool is_valid = (maxshare_map.find(pk) != maxshare_map.end());
                    send_out(serverfd, &is_valid, sizeof(bool));
                    if (!is_valid)
                        continue;
                    memcpy(share, maxshare_map[pk], share_sz);
                    for (int j = 0; j <= B; j++)
                        a[j] ^= share[j];
                }
                uint32_t b[B+1];
                for (int j = 0; j <= B; j++)
                    read_in(serverfd, &(b[j]), sizeof(uint32_t));
                for (int j = B; j >= 0; j--) {
                    // std::cout << a[j] << " " << b[j] << std::endl;
                    a[j] ^= b[j];
                    if (a[j] != 0) {
                        std::cout << "Maximum: " << j << std::endl;
                        break;
                    }
                }
            }
            long long t = time_from(start);
            std::cout << "Time taken : " << t << std::endl;
        } else if (msg.type == VAR_OP) {
            // Alternate idea: make each of these a function. Var calls intsum twice, then snip. 
            std::cout << "VAR_OP" << std::endl;
            auto start = clock_start();  // for benchmark

            std::vector<VarShare> varshares;
            std::unordered_map<std::string, uint32_t> varshare_map;
            std::unordered_map<std::string, uint32_t> varshare_map_squared;
            std::unordered_map<std::string, ClientPacket> varshare_map_packets;

            VarShare varshare;
            const uint32_t small_max = 1 << (num_bits / 2);  // for values
            const uint32_t var_max = 1 << num_bits;      // for values squared
            const size_t num_inputs = msg.num_of_inputs;

            std::cout << "small_max: " << small_max << std::endl;
            std::cout << "var_max: " << var_max << std::endl;

            ShareReceiver client_share_receiver(newsockfd);

            for (int i = 0; i < num_inputs; i++) {
                std::cout << "i = " << i << std::endl;
                int bytes_read = read_in(newsockfd, &varshare, sizeof(VarShare));
                if (bytes_read == 0) error_exit("Read 0 bytes. Connection closed?");
                std::string pk(varshare.pk, varshare.pk+32);

                ClientPacket packet = nullptr;
                client_share_receiver.client_packet(packet);

                if ((varshare_map.find(pk) != varshare_map.end())
                    or (varshare.val >= small_max) 
                    or (varshare.val_squared >= var_max)) {
                    std::cout << " invalid" << std::endl;
                    continue;  // Reject the input
                }

                std::cout << " valid" << std::endl;
                varshares.push_back(varshare);
                varshare_map[pk] = varshare.val;
                varshare_map_squared[pk] = varshare.val_squared;
                varshare_map_packets[pk] = packet;
            }

            std::cout << "Received " << num_inputs << " shares" << std::endl;
            std::cout << varshares.size() << " are valid" << std::endl;

            // Client shares complete. Communicate with each other.

            std::cout << "Forking for server chats" << std::endl;
            pid_t pid = fork();
            if (pid == 0) {
                std::cout << " Forked Child leaving." << std::endl;
                continue;
            }
            std::cout << " Parent for server chat" << std::endl;

            if (server_num == 1) {
                std::cout << "Server 1 chatting" << std::endl;
                size_t num_inputs = varshares.size();
                uint32_t shares[num_inputs];
                uint32_t shares_squared[num_inputs];

                std::cout << "sending num inputs: " << num_inputs << std::endl;
                // communicates number of valid on this side.
                send_out(serverfd, &num_inputs, sizeof(size_t));

                // Only run roots_init once on N.
                bool have_roots_init = false;

                for (int i = 0; i < num_inputs; i++) {
                    std::cout << "i = " << i << std::endl;
                    send_out(serverfd, &varshares[i].pk, 32);

                    bool other_valid;
                    read_in(serverfd, &other_valid, sizeof(bool));

                    std::cout << " got other valid: " << std::boolalpha << other_valid << std::endl;
                    if (!other_valid)
                        continue;

                    std::string pk(varshares[i].pk, varshares[i].pk+32);
                    shares[i] = varshares[i].val;
                    shares_squared[i] = varshares[i].val_squared;
                    ClientPacket packet = varshare_map_packets[pk];

                    // SNIPS
                    Circuit* circuit = CheckVar();
                    if (not have_roots_init) {
                        int N = NextPowerofTwo(circuit->NumMulGates() + 1);
                        init_roots(N);
                        have_roots_init = true;
                    }

                    bool circuit_valid = run_snip(circuit, packet, randomX, server_share_sender, server_share_receiver, server_num);
                    std::cout << " Circuit for " << i << " validity: " << std::boolalpha << circuit_valid << std::endl;
                    send_out(serverfd, &circuit_valid, sizeof(circuit_valid));
                }

                // Compute result
                NetIO* io = new NetIO(SERVER0_IP, 60051);
                uint64_t b = intsum_ot_receiver(io, &shares[0], num_inputs, num_bits);
                std::cout << "Sending b = " << b << std::endl;
                send_out(serverfd, &b, sizeof(uint64_t));

                uint64_t b2 = intsum_ot_receiver(io, &shares_squared[0], num_inputs, num_bits);
                std::cout << "Sending b2 = " << b << std::endl;
                send_out(serverfd, &b2, sizeof(uint64_t));
            } else {
                size_t num_inputs;
                std::cout << "waiting for num inputs" << std::endl;
                read_in(serverfd, &num_inputs, sizeof(size_t));

                std::cout << "num inputs: " << num_inputs << std::endl;

                uint32_t shares[num_inputs];
                uint32_t shares_squared[num_inputs];
                bool valid[num_inputs];

                bool have_roots_init = false;

                for (int i = 0; i < num_inputs; i++) {
                    std::cout << "i = " << i << std::endl;
                    char pk_buf[32];
                    read_in(serverfd, &pk_buf[0], 32);
                    std::string pk(pk_buf, pk_buf + 32);

                    bool is_valid = (varshare_map.find(pk) != varshare_map.end());
                    std::cout << " is_valid = " << is_valid << std::endl;

                    valid[i] = is_valid;
                    if (!is_valid)
                        continue;

                    // Communicate validity, to run snips less.
                    send_out(serverfd, &is_valid, sizeof(bool));

                    shares[i] = varshare_map[pk];
                    shares_squared[i] = varshare_map_squared[pk];
                    ClientPacket packet = varshare_map_packets[pk];

                    Circuit* circuit = CheckVar();
                    if (not have_roots_init) {
                        int N = NextPowerofTwo(circuit->NumMulGates() + 1);
                        init_roots(N);
                        have_roots_init = true;
                    }

                    bool circuit_valid = run_snip(circuit, packet, randomX, server_share_sender, server_share_receiver, server_num);
                    std::cout << " Circuit for " << i << " validity: " << std::boolalpha << circuit_valid << std::endl;
                    if (!circuit_valid) {
                        valid[i] = false;
                    }

                    bool other_valid;
                    read_in(serverfd, &other_valid, sizeof(bool));
                    std::cout << " Other circuit valid: " << std::boolalpha << other_valid << std::endl;
                    if (!other_valid) {
                        valid[i] = false;
                    }
                }

                // Compute result
                NetIO* io = new NetIO(nullptr, 60051);
                uint64_t a = intsum_ot_sender(io, &shares[0], &valid[0], num_inputs, num_bits);
                std::cout << "have a: " << a << std::endl;
                uint64_t a2 = intsum_ot_sender(io, &shares_squared[0], &valid[0], num_inputs, num_bits);
                std::cout << "have a2: " << a2 << std::endl;
                uint64_t b, b2;

                read_in(serverfd, &b, sizeof(uint64_t));
                std::cout << "got b: " << b << std::endl;
                read_in(serverfd, &b2, sizeof(uint64_t));
                std::cout << "got b2: " << b2 << std::endl;
                float ex = 1.0 * (a + b) / num_inputs;
                float ex2 = 1.0 * (a2 + b2) / num_inputs;
                float ans = ex2 - (ex * ex);
                std::cout << "Ans: " << ex2 << " - (" << ex << ")^2 = " << ans << std::endl;
            }
            long long t = time_from(start);
            std::cout << "Time taken : " << t << std::endl;
        } else {
            std::cout << "Unrecognized message type: " << msg.type << std::endl;
        }

        std::cout << "end of loop" << std::endl << std::endl;
        close(newsockfd);
    }

    return 0;
}
