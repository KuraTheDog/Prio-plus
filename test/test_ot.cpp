#include <iostream>

#include "utils_test_connect.h"
#include "../constants.h"
#include "../net_share.h"
#include "../ot.h"

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

  OT_Wrapper* ot0 = new OT_Wrapper(SERVER0_IP, SERVER0_OT_PORT, server_num == 0);
  OT_Wrapper* ot1 = new OT_Wrapper(SERVER1_IP, SERVER1_OT_PORT, server_num == 1);

  std::cout << "Simple ot test" << std::endl;

  if (server_num == 0) {
    uint64_t a[1] = {10};
    uint64_t b[1] = {20};
    ot0->send(a, b, 1);
  } else {
    uint64_t* d = new uint64_t[1];
    bool c[1] = {1};
    ot0->recv(d, c, 1);
  }

  std::cout << "Making bool triples" << std::endl;

  auto triples = gen_boolean_beaver_triples(server_num, m, ot0, ot1);

  if (pid == 0) {
    sleep(1);
    int cli_sockfd = init_sender();

    std::cout << "Validating bool triples" << std::endl;
    BooleanBeaverTriple* other_triple = new BooleanBeaverTriple();
    for (int i = 0; i < m; i++) {
      recv_BooleanBeaverTriple(cli_sockfd, other_triple);
      BooleanBeaverTriple* triple = triples.front();
      triples.pop();
      std::cout << "ab = (" << triple->a << " ^ " << other_triple->a << ") * (";
      std::cout << triple->b << " ^ " << other_triple->b << ") = ";
      std::cout << ((triple->a ^ other_triple->a) & (triple->b ^ other_triple->b));
      std::cout << ", vs c = " << triple->c << " ^ " << other_triple->c << " = " << (triple->c ^ other_triple->c) << std::endl;
    }
    delete other_triple;

    std::cout << "Making arith triple" << std::endl;
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

    delete other_btriple;
    delete btriple;

    close(cli_sockfd);
  } else {
    int sockfd = init_receiver();
    int newsockfd = accept_receiver(sockfd);

    for (int i = 0; i < m; i++) {
      BooleanBeaverTriple* triple = triples.front();
      triples.pop();
      send_BooleanBeaverTriple(newsockfd, triple);
      delete triple;
    }

    BeaverTriple* btriple = generate_beaver_triple_lazy(newsockfd, server_num);

    send_BeaverTriple(newsockfd, btriple);

    delete btriple;

    close(sockfd);
    close(newsockfd);
  }

  delete ot0;
  delete ot1;

  clear_constants();

  return 0;
}