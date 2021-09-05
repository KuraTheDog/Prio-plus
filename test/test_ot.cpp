#include <iostream>

#include "utils_test_connect.h"
#include "../constants.h"
#include "../net_share.h"
#include "../ot.h"

#define SERVER0_IP "127.0.0.1"
#define SERVER1_IP "127.0.0.1"

#define MOD 100
// #define MOD 0  // 2^64, default

void run_server0(const size_t m) {
  const int server_num = 0;

  int cli_sockfd = init_sender();

  OT_Wrapper* ot0 = new OT_Wrapper(server_num == 0 ? nullptr : SERVER0_IP, 60051);
  OT_Wrapper* ot1 = new OT_Wrapper(server_num == 1 ? nullptr : SERVER1_IP, 60052);

  // gen beaver triple
  // auto triples = gen_boolean_beaver_triples(server_num, m, ot0, ot1);

  // std::cout << "Validating bool triples" << std::endl;
  // BooleanBeaverTriple* other_triple = new BooleanBeaverTriple();
  // for (int i = 0; i < m; i++) {
  //   recv_BooleanBeaverTriple(cli_sockfd, other_triple);
  //   BooleanBeaverTriple* triple = triples.front();
  //   triples.pop();
  //   std::cout << "ab = (" << triple->a << " ^ " << other_triple->a << ") * (";
  //   std::cout << triple->b << " ^ " << other_triple->b << ") = ";
  //   std::cout << ((triple->a ^ other_triple->a) & (triple->b ^ other_triple->b));
  //   std::cout << ", vs c = " << triple->c << " ^ " << other_triple->c << " = " << (triple->c ^ other_triple->c) << std::endl;
  // }
  // delete other_triple;

  // // gen arith triple

  // std::cout << "Making arith triple" << std::endl;
  // BeaverTriple* btriple = generate_beaver_triple_lazy(cli_sockfd, server_num);
  // std::cout << "Validating lazy arith triple" << std::endl;
  // BeaverTriple* other_btriple = new BeaverTriple();
  // recv_BeaverTriple(cli_sockfd, other_btriple);

  // fmpz_t tmp; fmpz_init(tmp);
  // fmpz_t tmp2; fmpz_init(tmp2);
  // fmpz_add(tmp, btriple->A, other_btriple->A);
  // fmpz_mod(tmp, tmp, Int_Modulus);
  // fmpz_add(tmp2, btriple->B, other_btriple->B);
  // fmpz_mod(tmp2, tmp2, Int_Modulus);
  // fmpz_mul(tmp, tmp, tmp2);
  // fmpz_mod(tmp, tmp, Int_Modulus);
  // std::cout << "actual product: ("; fmpz_print(btriple->A);
  // std::cout << " + "; fmpz_print(other_btriple->A);
  // std::cout << ") * ("; fmpz_print(btriple->B);
  // std::cout << " + "; fmpz_print(other_btriple->B);
  // std::cout << ") = "; fmpz_print(tmp);
  // std::cout << std::endl;
  // fmpz_add(tmp, btriple->C, other_btriple->C);
  // fmpz_mod(tmp, tmp, Int_Modulus);
  // fmpz_print(btriple->C); std::cout << " + "; fmpz_print(other_btriple->C);
  // std::cout << " = "; fmpz_print(tmp); std::cout << std::endl;

  // delete other_btriple;
  // delete btriple;

  // Bitsum OT
  const bool valid[1] = {true};

  std::cout << "Testing bit OT convert" << std::endl;
  const bool bitshares[1] = {(m / 2) % 2 ? true : false};
  std::cout << "bit share0[0] = " << bitshares[0] << std::endl;
  const uint64_t a = bitsum_ot_sender(ot0, bitshares, valid, 1, MOD);
  uint64_t b;
  recv_uint64(cli_sockfd, b);
  std::cout << "bit have a = " << a << std::endl;
  std::cout << "bit got  b = " << b << std::endl;
  std::cout << "bit ans: " << a + b << std::endl;

  // Intsum OT

  std::cout << "Testing int OT convert" << std::endl;
  // const uint64_t intshares[1][1] = {{m ^ 9}};
  uint64_t** const intshares = new uint64_t*[1];
  intshares[0] = new uint64_t[1];
  intshares[0][0] = m ^ 9;
  const size_t sizes[1] = {4};
  std::cout << "int share0[0] = " << intshares[0][0] << std::endl;
  uint64_t** int_a = intsum_ot_sender(ot0, intshares, valid, sizes, 1, 1, MOD);
  delete[] intshares[0];
  delete[] intshares;
  recv_uint64(cli_sockfd, b);
  std::cout << "int have a = " << int_a[0][0] << std::endl;
  std::cout << "int got  b = " << b << std::endl;
  std::cout << "int ans: " << int_a[0][0] + b << std::endl;
  delete[] int_a[0];
  delete[] int_a;

  // cleanup

  delete ot0;
  delete ot1;
  close(cli_sockfd);
}

void run_server1(const size_t m) {
  const int server_num = 1;

  int sockfd = init_receiver();
  int newsockfd = accept_receiver(sockfd);

  OT_Wrapper* ot0 = new OT_Wrapper(server_num == 0 ? nullptr : SERVER0_IP, 60051);
  OT_Wrapper* ot1 = new OT_Wrapper(server_num == 1 ? nullptr : SERVER1_IP, 60052);

  // auto triples = gen_boolean_beaver_triples(server_num, m, ot0, ot1);
  // for (int i = 0; i < m; i++) {
  //   BooleanBeaverTriple* triple = triples.front();
  //   triples.pop();
  //   send_BooleanBeaverTriple(newsockfd, triple);
  //   delete triple;
  // }

  // BeaverTriple* btriple = generate_beaver_triple_lazy(newsockfd, server_num);
  // send_BeaverTriple(newsockfd, btriple);
  // delete btriple;

  const bool bitshares[1] = {m % 2 ? true : false};
  std::cout << "bit share1[0] = " << bitshares[0] << std::endl;
  const uint64_t b = bitsum_ot_receiver(ot0, bitshares, 1, MOD);
  send_uint64(newsockfd, b);

  uint64_t** const intshares = new uint64_t*[9];
  intshares[0] = new uint64_t[1];
  intshares[0][0] = 9;
  const size_t sizes[1] = {4};
  std::cout << "int share1[0] = " << intshares[0][0] << std::endl;
  uint64_t** int_b = intsum_ot_receiver(ot0, intshares, sizes, 1, 1, MOD);
  delete[] intshares[0];
  delete[] intshares;
  send_uint64(newsockfd, int_b[0][0]);
  delete[] int_b[0];
  delete[] int_b;

  delete ot0;
  delete ot1;
  close(sockfd);
  close(newsockfd);
}

int main(int argc, char** argv){
  int m = 4;
  if(argc >= 2){
    m = atoi(argv[1]);
  }

  init_constants();

  pid_t pid = fork();

  if (pid == 0) {
    run_server0(m);
  } else {
    run_server1(m);
  }

  clear_constants();

  return 0;
}