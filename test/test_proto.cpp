#include <iostream>

#include "utils_test_connect.h"
#include "../constants.h"
#include "../net_share.h"
#include "../proto.h"

#define SERVER0_IP "127.0.0.1"
#define SERVER1_IP "127.0.0.1"

int main(int argc, char** argv){
  int m = 10;
  if(argc >= 2){
    m = atoi(argv[1]);
  }

  init_constants();

  pid_t pid = fork();
  int server_num = (pid == 0 ? 0 : 1);

  std::cout << "Making io objects" << std::endl;

  NetIO* io0 = new NetIO(server_num == 0 ? nullptr : SERVER0_IP, 60051, true);
  NetIO* io1 = new NetIO(server_num == 1 ? nullptr : SERVER1_IP, 60052, true); 

  std::cout << "Making bool triples" << std::endl;

  BooleanBeaverTriple* triples = gen_boolean_beaver_triples(server_num, m, io0, io1);

  if (pid == 0) {
    sleep(1);
    int cli_sockfd = init_sender();

    std::cout << "Validating bool triples" << std::endl;
    BooleanBeaverTriple triple;
    for (int i = 0; i < m; i++) {
      recv_BooleanBeaverTriple(cli_sockfd, triple);
      std::cout << "ab = (" << triple.a << " ^ " << triples[i].a << ") * (";
      std::cout << triple.b << " ^ " << triples[i].b << ") = ";
      std::cout << ((triple.a ^ triples[i].a) & (triple.b ^ triples[i].b));
      std::cout << ", vs c = " << triple.c << " ^ " << triples[i].c << " = " << (triple.c ^ triples[i].c) << std::endl;
    }

    std::cout << "Making arith triple" << std::endl;
    // BeaverTriple* btriple = generate_beaver_triple(cli_sockfd, server_num, io0, io1);
    BeaverTriple* btriple = generate_beaver_triple_lazy(cli_sockfd, server_num);
    std::cout << "Validating arith triple" << std::endl;
    BeaverTriple* other_btriple = new BeaverTriple();
    recv_BeaverTriple(cli_sockfd, other_btriple);

    fmpz_t tmp; fmpz_init(tmp);
    fmpz_t tmp2; fmpz_init(tmp2);
    fmpz_add(tmp, btriple->A, other_btriple->A);
    fmpz_mod(tmp, tmp, Int_Modulus);
    fmpz_add(tmp2, btriple->B, other_btriple->B);
    fmpz_mod(tmp2, tmp2, Int_Modulus);
    fmpz_mul(tmp, tmp, tmp2);
    fmpz_mod(tmp, tmp, Int_Modulus);
    std::cout << "actual product: ("; fmpz_print(btriple->A);
    std::cout << " + "; fmpz_print(other_btriple->A);
    std::cout << ") * ("; fmpz_print(btriple->B);
    std::cout << " + "; fmpz_print(other_btriple->B);
    std::cout << ") = "; fmpz_print(tmp);
    std::cout << std::endl;
    fmpz_add(tmp, btriple->C, other_btriple->C);
    fmpz_mod(tmp, tmp, Int_Modulus);
    fmpz_print(btriple->C); std::cout << " + "; fmpz_print(other_btriple->C);
    std::cout << " = "; fmpz_print(tmp); std::cout << std::endl;

    close(cli_sockfd);
  } else {
    int sockfd = init_receiver();
    int newsockfd = accept_receiver(sockfd);

    for (int i = 0; i < m; i++)
      send_BooleanBeaverTriple(newsockfd, triples[i]);

    // BeaverTriple* btriple = generate_beaver_triple(newsockfd, server_num, io0, io1);
    BeaverTriple* btriple = generate_beaver_triple_lazy(newsockfd, server_num);

    send_BeaverTriple(newsockfd, btriple);

    close(sockfd);
    close(newsockfd);
  }

  delete io0;
  delete io1;

  return 0;
}