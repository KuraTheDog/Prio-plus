#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <unistd.h>

#include <iostream>

#include "../constants.h"
#include "../net_share.h"
#include "../server.h"
#include "../share.h"

int PORT = 8888;

void error_exit(const char *msg){
  perror(msg);
  exit(EXIT_FAILURE);
}

int init_sender() {
  std::cout << "send: start" << std::endl;

  int sockfd;
  struct sockaddr_in addr;

  sockfd = socket(AF_INET, SOCK_STREAM, 0);
  if (sockfd < 0) error_exit("send: ERROR opening socket");

  int reuse = 1;
  if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, &reuse, sizeof(reuse)))
    error_exit("Sockopt failed");
  if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEPORT, &reuse, sizeof(reuse)))
    error_exit("Sockopt failed");

  bzero((char *) &addr, sizeof(addr));
  addr.sin_port = htons(PORT);
  addr.sin_family = AF_INET;
  inet_pton(AF_INET, "127.0.0.1", &addr.sin_addr);

  if (connect(sockfd, (sockaddr*) &addr, sizeof(addr)) < 0)
      error_exit("ERROR on connect");

  return sockfd;
}

int init_receiver() {
  std::cout << "recv: start" << std::endl;

  int sockfd;
  struct sockaddr_in rec_addr;

  sockfd = socket(AF_INET, SOCK_STREAM, 0);
  if (sockfd < 0) error_exit("recv: ERROR opening socket");

  int reuse = 1;
  if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, &reuse, sizeof(reuse)))
    error_exit("Sockopt failed");
  if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEPORT, &reuse, sizeof(reuse)))
    error_exit("Sockopt failed");

  bzero((char *) &rec_addr, sizeof(rec_addr));
  rec_addr.sin_family = AF_INET;
  rec_addr.sin_addr.s_addr = INADDR_ANY;
  rec_addr.sin_port = htons(PORT);

  if (bind(sockfd, (struct sockaddr *) &rec_addr, sizeof(rec_addr)) < 0)
    error_exit("recv: ERROR on binding");

  if(listen(sockfd, 2) < 0) error_exit("recv: ERROR on listen");
  return sockfd;
}

int accept_receiver(int sockfd) {
  int newsockfd;
  socklen_t snd_len;
  struct sockaddr_in snd_addr;
  snd_len = sizeof(snd_addr);
  newsockfd = accept(sockfd, (struct sockaddr *) & snd_addr, &snd_len);
  if (newsockfd < 0) error_exit("ERROR on accept");

  printf("recv: got connection from %s port %d\n",
         inet_ntoa(snd_addr.sin_addr), ntohs(snd_addr.sin_port));
  return newsockfd;
}

void localTest(const size_t n) {
  std::cout << "Running Local Test" << std::endl;

  fmpz_t tmp;
  fmpz_init(tmp);

  // random adjusting. comment to change seed.
  // fmpz_randm(tmp, seed, Int_Modulus);
  // fmpz_randm(tmp, seed, Int_Modulus);

  EdaBit* ebit00 = new EdaBit(n);  // 0's share of 0
  EdaBit* ebit01 = new EdaBit(n);  // 1's share of 0
  makeLocalEdaBit(ebit00, ebit01, n);

  EdaBit* ebit10 = new EdaBit(n);  // 0's share of 1
  EdaBit* ebit11 = new EdaBit(n);  // 1's share of 1
  makeLocalEdaBit(ebit10, ebit11, n);

  // 0 has ebit00, ebit01. 1 has ebit10, ebit11
  // They move things around
  // 0 has ebit00, ebit10. 1 has ebit01, ebit11

  // std::cout << "0 has ebit00: \n"; ebit00->print();
  // std::cout << "0 has ebit10: \n"; ebit10->print();

  // std::cout << "1 has ebit01: \n"; ebit01->print();
  // std::cout << "1 has ebit11: \n"; ebit11->print();

  // value checking
  // fmpz_add(tmp, ebit00->r, ebit01->r);
  // fmpz_mod(tmp, tmp, Int_Modulus);
  // std::cout << "ebit00 and ebit01, sum = "; fmpz_print(tmp); std::cout << std::endl;
  // std::cout << "ebit00 and ebit01, xor = " << (ebit00->get_int_b() ^ ebit01->get_int_b()) << std::endl;

  // fmpz_add(tmp, ebit10->r, ebit11->r);
  // fmpz_mod(tmp, tmp, Int_Modulus);
  // std::cout << "ebit10 and ebit11, sum = "; fmpz_print(tmp); std::cout << std::endl;
  // std::cout << "ebit10 and ebit11, xor = " << (ebit10->get_int_b() ^ ebit11->get_int_b()) << std::endl;

  // Final bits
  EdaBit* ebit0 = new EdaBit(n);
  EdaBit* ebit1 = new EdaBit(n);

  // r'
  fmpz_add(ebit0->r, ebit00->r, ebit10->r);
  fmpz_mod(ebit0->r, ebit0->r, Int_Modulus);
  // std::cout << "[r']^0 = "; fmpz_print(ebit0->r); std::cout << std::endl;

  fmpz_add(ebit1->r, ebit01->r, ebit11->r);
  fmpz_mod(ebit1->r, ebit1->r, Int_Modulus);
  // std::cout << "[r']^1 = "; fmpz_print(ebit1->r); std::cout << std::endl;

  // bs
  // x0 = ebit00->b, y0 = ebit10->b
  // x1 = ebit00->b, y1 = ebit10->b
  // z0 = ebit0->b, z1 = ebit1->b
  // c_{i+1} = c_i xor ((x_i xor c_i) and (y_i xor c_i))
  // output z_i = x_i xor y_i xor c_i

  bool carry0 = false, carry1 = false;  // carry bits
  bool x0, x1, y0, y1, tmp0, tmp1;  // addition calculation
  bool d, e, d0, d1, e0, e1;  // triple calculation
  bool a0, a1, b0, b1, c0, c1;  // beaver triple;

  for (int i = 0; i < n; i++) {
    ebit0->b[i] = carry0 ^ ebit00->b[i] ^ ebit10->b[i];
    x0 = carry0 ^ ebit00->b[i];
    y0 = carry0 ^ ebit10->b[i];

    ebit1->b[i] = carry1 ^ ebit01->b[i] ^ ebit11->b[i];
    x1 = carry1 ^ ebit01->b[i];
    y1 = carry1 ^ ebit11->b[i];

    // boolean beaver triple.
    // (a0^a1) * (b0^b1) = (c0^c1). 5 random per.
    fmpz_randbits(tmp, seed, 5);
    a0 = get_fmpz_bit(tmp, 0);
    a1 = get_fmpz_bit(tmp, 1);
    b0 = get_fmpz_bit(tmp, 2);
    b1 = get_fmpz_bit(tmp, 3);
    c0 = get_fmpz_bit(tmp, 4);
    c1 = c0 ^ ((a0 ^ a1) and (b0 ^ b1));

    // [c] xor [x]*e xor [y] * d xor e * d
    
    d0 = x0 ^ a0;
    e0 = y0 ^ b0;

    d1 = x1 ^ a1;
    e1 = y1 ^ b1;
    // broadcast
    d = d0 ^ d1;
    e = e0 ^ e1;

    tmp0 = c0 ^ (x0 and e) ^ (y0 and d) ^ (d and e);
    carry0 ^= tmp0;

    tmp1 = c1 ^ (x1 and e) ^ (y1 and d);
    carry1 ^= tmp1;
  }

  // std::cout << " [b]_2^0 = " << ebit0->get_int_b() << (carry0? " with carry" : " no carry") << std::endl;
  // std::cout << " [b]_2^1 = " << ebit1->get_int_b() << (carry1? " with carry" : " no carry") << std::endl;

  // b_p's

  // dabit for adjusting. Doing one bit convertToField
  // TODO: actually do via servers
  DaBit* dabit0 = new DaBit();
  DaBit* dabit1 = new DaBit();
  makeLocalDaBit(dabit0, dabit1);
  // std::cout << "conversion dabit0\n"; dabit0->print();
  // std::cout << "conversion dabit1\n"; dabit1->print();
  bool v0 = carry0 ^ dabit0->b2;
  bool v1 = carry1 ^ dabit1->b2;
  // 0 has v0, 1 has v1. Broadcast both so both know v.
  bool v = v0 ^ v1;
  fmpz_t bp0;
  if (v) {  // v = 1, so x0 = 1 - [b]_p
    fmpz_init_set_si(bp0, 1);
    fmpz_sub(bp0, bp0, dabit0->bp);
    fmpz_mod(bp0, bp0, Int_Modulus);
  } else {  // v = 0, so x0 = [b]_p
    fmpz_init_set(bp0, dabit0->bp);
  }
  // std::cout << "bp0 = "; fmpz_print(bp0); std::cout << std::endl;
  fmpz_t bp1;
  if (v) {  // v = 1, so x1 = [b]_p
    fmpz_init_set_si(bp1, 0);
    fmpz_sub(bp1, bp1, dabit1->bp);
    fmpz_mod(bp1, bp1, Int_Modulus);
  } else {  // v = 0, so x1 = [b]_p
    fmpz_init_set(bp1, dabit1->bp);
  }
  // std::cout << "bp1 = "; fmpz_print(bp1); std::cout << std::endl;

  // r = r' - 2^n * bp
  fmpz_mul_ui(tmp, bp0, 1 << n);
  fmpz_sub(ebit0->r, ebit0->r, tmp);
  fmpz_mod(ebit0->r, ebit0->r, Int_Modulus);

  fmpz_mul_ui(tmp, bp1, 1 << n);
  fmpz_sub(ebit1->r, ebit1->r, tmp);
  fmpz_mod(ebit1->r, ebit1->r, Int_Modulus);

  std::cout << "Final ebit0:\n"; ebit0->print();
  std::cout << "Final ebit1:\n"; ebit1->print();

  fmpz_add(tmp, ebit0->r, ebit1->r);
  fmpz_mod(tmp, tmp, Int_Modulus);
  std::cout << "Validation" << std::endl;
  std::cout << " sum = "; fmpz_print(tmp); std::cout << std::endl;
  std::cout << " xor = " << (ebit0->get_int_b() ^ ebit1->get_int_b()) << std::endl;
}

void runServerTest(const int server_num, const int serverfd, const size_t n) {

  // TODO: use proper triple sharing
  BooleanBeaverTriple triples[n];
  for (int i = 0; i < n; i++) {
    if (server_num == 0) {
      fmpz_t tmp;
      fmpz_init(tmp);
      fmpz_randbits(tmp, seed, 5);
      triples[i].a = get_fmpz_bit(tmp, 0);
      triples[i].b = get_fmpz_bit(tmp, 1);
      triples[i].c = get_fmpz_bit(tmp, 2);
      bool a2 = get_fmpz_bit(tmp, 3);
      bool b2 = get_fmpz_bit(tmp, 4);
      bool c2 = triples[i].c ^ ((triples[i].a ^ a2) and (triples[i].b ^ b2));

      send_bool(serverfd, a2);
      send_bool(serverfd, b2);
      send_bool(serverfd, c2);
      fmpz_clear(tmp);
    } else {
      recv_bool(serverfd, triples[i].a);
      recv_bool(serverfd, triples[i].b);
      recv_bool(serverfd, triples[i].c);
    }
  }

  // TODO: Use proper dabit generation
  DaBit* dabit = new DaBit();
  if (server_num == 0) {
    DaBit* dabit_other = new DaBit();
    makeLocalDaBit(dabit, dabit_other);
    send_DaBit(serverfd, dabit_other);
    delete dabit_other;
  } else {
    recv_DaBit(serverfd, dabit);
  }

  EdaBit* ebit = generateEdaBit(serverfd, server_num, n, triples, dabit);

  std::cout << "Final ebit" << server_num << ":\n";
  ebit->print();

  // Validation
  std::cout << "Validation" << std::endl;
  if (server_num == 0) {
    send_EdaBit(serverfd, ebit, n);
  } else {
    EdaBit* ebit_other = new EdaBit(n);
    recv_EdaBit(serverfd, ebit_other, n);
    fmpz_t tmp;
    fmpz_init(tmp);

    fmpz_add(tmp, ebit->r, ebit_other->r);
    fmpz_mod(tmp, tmp, Int_Modulus);
    std::cout << " sum = "; fmpz_print(tmp); std::cout << std::endl;
    std::cout << " xor = " << (ebit->get_int_b() ^ ebit_other->get_int_b()) << std::endl;

    fmpz_clear(tmp);
    delete ebit_other;
  }

  delete dabit;
}

void serverTest(const size_t n) {
  std::cout << "Running server test" << std::endl;
  int sockfd = init_receiver();

  pid_t pid = fork();
  if (pid == 0) {
    int cli_sockfd = init_sender();

    fmpz_t tmp;
    fmpz_init(tmp);
    // fmpz_randm(tmp, seed, Int_Modulus);
    // fmpz_randm(tmp, seed, Int_Modulus);
    fmpz_clear(tmp);

    runServerTest(0, cli_sockfd, n);
    close(cli_sockfd);
  } else if (pid > 0) {
    int newsockfd = accept_receiver(sockfd);

    // alter randomness to be different from the sender
    fmpz_t tmp;
    fmpz_init(tmp);
    fmpz_randm(tmp, seed, Int_Modulus);
    fmpz_clear(tmp);

    runServerTest(1, newsockfd, n);

    close(newsockfd);
    close(sockfd);
  }
}

int main(int argc, char* argv[]) {
  init_constants();
  
  // Should be < log_2 Int_Modulus
  size_t n = 16;
  std::cout << "n = " << n << std::endl;

  localTest(n);
  std::cout << std::endl;
  serverTest(n);

  clear_constants();
  return 0;
}