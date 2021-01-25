#include <iostream>

#include "utils_test_connect.h"
#include "../net_share.h"
#include "../proto.h"

#define SERVER0_IP "127.0.0.1"
#define SERVER1_IP "127.0.0.1"

int main(int argc, char** argv){
  int m = 10;
  if(argc >= 2){
    m = atoi(argv[1]);
  }

  BooleanBeaverTriple* triples;

  pid_t pid = fork();
  int server_num = (pid == 0 ? 0 : 1);
  triples = gen_boolean_beaver_triples(server_num, m);  

  int sockfd = init_receiver();
  if (pid == 0) {
    int cli_sockfd = init_sender();

    BooleanBeaverTriple triple;
    for (int i = 0; i < m; i++) {
      recv_BooleanBeaverTriple(cli_sockfd, triple);
      std::cout << "ab = (" << triple.a << " ^ " << triples[i].a << ") * (";
      std::cout << triple.b << " ^ " << triples[i].b << ") = ";
      std::cout << ((triple.a ^ triples[i].a) & (triple.b ^ triples[i].b));
      std::cout << ", vs c = " << triple.c << " ^ " << triples[i].c << " = " << (triple.c ^ triples[i].c) << std::endl;
    }

    close(cli_sockfd);
  } else {
    int newsockfd = accept_receiver(sockfd);

    for (int i = 0; i < m; i++)
      send_BooleanBeaverTriple(newsockfd, triples[i]);

    close(newsockfd);
  }
  close(sockfd);

  return 0;
}