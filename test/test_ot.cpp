#include <iostream>

#include "utils_test_connect.h"
#include "../constants.h"
#include "../net_share.h"
#include "../ot.h"

#define SERVER0_IP "127.0.0.1"
#define SERVER1_IP "127.0.0.1"

int main(int argc, char** argv){
  int n = 10;
  if(argc >= 2){
    n = atoi(argv[1]);
  }
  int m = 1;
  if (argc >= 3)
    m = atoi(argv[2]);

  init_constants();

  pid_t pid = fork();
  int server_num = (pid == 0 ? 0 : 1);


  // TODO: serverfd stuff
  int serverfd = -1;

  #if OT_TYPE == LIBOTE_SILENT
  std::cout << "setting up sockets" << std::endl;
  if (server_num == 0) {
    int sockfd = init_receiver();
    int newsockfd = accept_receiver(sockfd);
    serverfd = newsockfd;
  } else if (server_num == 1) {
    sleep(1);
    int cli_sockfd = init_sender();
    serverfd = cli_sockfd;
  }
  #endif

  std::cout << "Making io objects" << std::endl;

  std::cout << "making ot0" << std::endl;
  OT_Wrapper* ot0 = new OT_Wrapper(SERVER0_IP, SERVER0_OT_PORT, server_num == 0, serverfd);
  std::cout << "making ot1" << std::endl;
  OT_Wrapper* ot1 = new OT_Wrapper(SERVER1_IP, SERVER1_OT_PORT, server_num == 1, serverfd);

  std::cout << "Simple ot test, n = " << n << ", iterations: " << m << std::endl;
  for (int j = 0; j < m; j++) {
    if (server_num == 0) {
      uint64_t a[n];
      uint64_t b[n];
      for (int i = 0; i < n; i++) {
        a[i] = (1000 * j + 10 * i + 1);
        b[i] = 10 * (1000 * j + 10 * i + 1);
        if (i < 4 or i == n-1)
          std::cout << j << ", " << i << " send " << a[i] << ", " << b[i] << std::endl;
      }
      ot0->send(a, b, n);
    } else {
      uint64_t* d = new uint64_t[n];
      bool c[n];
      for (int i = 0; i < n; i++)
        c[i] = i % 2;
      ot0->recv(d, c, n);
      for (int i = 0; i < n; i++) {
        if (i < 4 or i == n-1)
          std::cout << j << ", " << i << " got " << c[i] << " = " << d[i] << std::endl;
      }
    }
  }

  return 0;

  std::cout << "Making bool triples" << std::endl;

  auto triples = gen_boolean_beaver_triples(server_num, n, ot0, ot1);

  if (pid == 0) {
    sleep(1);
    int cli_sockfd = init_sender();

    std::cout << "Validating bool triples" << std::endl;
    BooleanBeaverTriple* other_triple = new BooleanBeaverTriple();
    for (int i = 0; i < n; i++) {
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

    for (int i = 0; i < n; i++) {
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