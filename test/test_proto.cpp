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

  BooleanBeaverTriple* triples = gen_boolean_beaver_triples(server_num, m);  

  int sockfd = init_receiver();
  if (pid == 0) {
    // random adjust
    fmpz_t jnk; fmpz_init(jnk);
    fmpz_randm(jnk, seed, Int_Modulus);

    int cli_sockfd = init_sender();

    BooleanBeaverTriple triple;
    for (int i = 0; i < m; i++) {
      recv_BooleanBeaverTriple(cli_sockfd, triple);
      std::cout << "ab = (" << triple.a << " ^ " << triples[i].a << ") * (";
      std::cout << triple.b << " ^ " << triples[i].b << ") = ";
      std::cout << ((triple.a ^ triples[i].a) & (triple.b ^ triples[i].b));
      std::cout << ", vs c = " << triple.c << " ^ " << triples[i].c << " = " << (triple.c ^ triples[i].c) << std::endl;
    }

    BeaverTriple* btriple = generate_beaver_triple(server_num, cli_sockfd);
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
    int newsockfd = accept_receiver(sockfd);

    for (int i = 0; i < m; i++)
      send_BooleanBeaverTriple(newsockfd, triples[i]);

    BeaverTriple* btriple = generate_beaver_triple(server_num, newsockfd);

    send_BeaverTriple(newsockfd, btriple);

    close(newsockfd);
  }
  close(sockfd);

  return 0;
}