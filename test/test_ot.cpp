#include <iostream>

#include "utils_test_connect.h"
#include "../constants.h"
#include "../net_share.h"
#include "../ot.h"

#define SERVER0_IP "127.0.0.1"
#define SERVER1_IP "127.0.0.1"

void run(const int server_num, const int sockfd, const int n, const int m, const char* addr0, const char* addr1) {

  std::cout << "Making io objects" << std::endl;

  auto start = clock_start();

  size_t batch_size = 1024;
  while((unsigned) n > batch_size)
    batch_size *= 2;

  OT_Wrapper* ot0 = new OT_Wrapper(addr0, SERVER0_OT_PORT, server_num == 0, sockfd, batch_size);
  std::cout << "Make ot0 time: " << sec_from(start) << std::endl;
  start = clock_start();
  OT_Wrapper* ot1 = new OT_Wrapper(addr1, SERVER1_OT_PORT, server_num == 1, sockfd, batch_size);
  std::cout << "Make ot1 time: " << sec_from(start) << std::endl;
  start = clock_start();

  #if OT_TYPE == LIBOTE_SILENT
  ot0->maybeUpdate();
  std::cout << "ot0 precompute time: " << sec_from(start) << std::endl;
  start = clock_start();
  ot1->maybeUpdate();
  std::cout << "ot1 precompute time: " << sec_from(start) << std::endl;
  start = clock_start();
  #endif

  std::cout << "Simple ot test, n = " << n << ", iterations: " << m << std::endl;
  for (int j = 0; j < m; j++) {
    if (server_num == 0) {
      uint64_t a[n];
      uint64_t b[n];
      for (int i = 0; i < n; i++) {
        a[i] = (10000 * j + 100 * i + 1);
        b[i] = 10 * (10000 * j + 100 * i + 1);
      }

      ot0->send(a, b, n);

      for (int i = 0; i < n; i++) {
        if (i < 3 or i == n-1)
          std::cout << j << ", " << i << " send " << a[i] << ", " << b[i] << std::endl;
      }
    } else {
      uint64_t* d = new uint64_t[n];
      bool c[n];
      for (int i = 0; i < n; i++)
        c[i] = i % 2;

      ot0->recv(d, c, n);

      for (int i = 0; i < n; i++) {
        if (i < 3 or i == n-1)
          std::cout << j << ", " << i << " got " << c[i] << " = " << d[i] << std::endl;
      }
    }
  }

  std::cout << "OT test time: " << sec_from(start) << std::endl;

  // sleep(100);
  return;

  std::cout << "Making bool triples" << std::endl;

  auto triples = gen_boolean_beaver_triples(server_num, n, ot0, ot1);

  if (server_num == 0) {
    std::cout << "Validating bool triples" << std::endl;
    BooleanBeaverTriple* other_triple = new BooleanBeaverTriple();
    for (int i = 0; i < n; i++) {
      recv_BooleanBeaverTriple(sockfd, other_triple);
      BooleanBeaverTriple* triple = triples.front();
      triples.pop();
      std::cout << "ab = (" << triple->a << " ^ " << other_triple->a << ") * (";
      std::cout << triple->b << " ^ " << other_triple->b << ") = ";
      std::cout << ((triple->a ^ other_triple->a) & (triple->b ^ other_triple->b));
      std::cout << ", vs c = " << triple->c << " ^ " << other_triple->c << " = " << (triple->c ^ other_triple->c) << std::endl;
    }
    delete other_triple;

    std::cout << "Making arith triple" << std::endl;
    BeaverTriple* btriple = generate_beaver_triple_lazy(sockfd, server_num);
    std::cout << "Validating arith triple" << std::endl;
    BeaverTriple* other_btriple = new BeaverTriple();
    recv_BeaverTriple(sockfd, other_btriple);

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
  } else {
    for (int i = 0; i < n; i++) {
      BooleanBeaverTriple* triple = triples.front();
      triples.pop();
      send_BooleanBeaverTriple(sockfd, triple);
      delete triple;
    }

    BeaverTriple* btriple = generate_beaver_triple_lazy(sockfd, server_num);

    send_BeaverTriple(sockfd, btriple);

    delete btriple;
  }

  delete ot0;
  delete ot1;
}


int main(int argc, char** argv) {

  int server_num = -1;
  if(argc >= 2)
    server_num = atoi(argv[1]);

  if (server_num < -1 or server_num > 1) {
    perror("Server num must be -1 (for fork), 0, or 1");
    exit(1);
  }

  // batch size
  int n = 10;
  if(argc >= 3)
    n = atoi(argv[2]);

  // num batches
  int m = 1;
  if (argc >= 4)
    m = atoi(argv[3]);

  init_constants();

  int sockfd = -1;
  if (server_num == -1 or server_num == 0) {
    sockfd = init_receiver();
  }

  const char* addr0 = SERVER0_IP;
  const char* addr1 = SERVER1_IP;

  if (server_num == -1) {
    addr0 = "127.0.0.1";
    addr1 = "127.0.0.1";
    pid_t pid = fork();
    server_num = (pid == 0 ? 0 : 1);
    if (server_num == 1)
      sleep(1);
  }

  if (server_num == 0) {
    int newsockfd = accept_receiver(sockfd);
    run(0, newsockfd, n, m, addr0, addr1);
    close(newsockfd);
    close(sockfd);
  } else {
    int cli_sockfd = init_sender(addr0);
    run(1, cli_sockfd, n, m, addr0, addr1);
    close(cli_sockfd);
  }

  clear_constants();
}
