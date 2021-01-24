#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <unistd.h>

#include <iostream>

#include "../constants.h"
#include "../edabit.h"
#include "../net_share.h"
#include "../proto.h"
#include "../share.h"

int PORT = 8887;

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

// Test version
void makeBeaverTripleLocally(fmpz_t a0, fmpz_t a1, fmpz_t b0, fmpz_t b1, fmpz_t c0, fmpz_t c1) {
  fmpz_init(a0);
  fmpz_randm(a0, seed, Int_Modulus);
  fmpz_init(a1);
  fmpz_randm(a1, seed, Int_Modulus);
  fmpz_init(b0);
  fmpz_randm(b0, seed, Int_Modulus);
  fmpz_init(b1);
  fmpz_randm(b1, seed, Int_Modulus);
  fmpz_init(c0);
  fmpz_randm(c0, seed, Int_Modulus);

  // c1 = (a0 + a1)(b0 + b1) - c0
  fmpz_t tmp;
  fmpz_init(c1);  
  fmpz_add(c1, a0, a1);
  fmpz_init(tmp);
  fmpz_add(tmp, b0, b1);
  fmpz_mul(c1, c1, tmp);
  fmpz_sub(c1, c1, c0);
  fmpz_mod(c1, c1, Int_Modulus);

  fmpz_clear(tmp);
}

// Local version of generateDaBit
void makeSharedDabitLocally(DaBit* dabit0, DaBit* dabit1) {
  DaBit* dabit00 = new DaBit();
  DaBit* dabit01 = new DaBit();
  makeLocalDaBit(dabit00, dabit01);
  DaBit* dabit10 = new DaBit();
  DaBit* dabit11 = new DaBit();
  makeLocalDaBit(dabit10, dabit11);

  // do xors. 
  // binary straightforward
  dabit0->b2 = dabit00->b2 ^ dabit10->b2;
  dabit1->b2 = dabit01->b2 ^ dabit11->b2;

  // a xor b = a + b - 2ab
  fmpz_t prod0, prod1;
  fmpz_init(prod0);
  fmpz_init(prod1);

  // (a0 + a1)(b0 + b1) = (c0 + c1)
  fmpz_t a0, a1, b0, b1, c0, c1;
  makeBeaverTripleLocally(a0, a1, b0, b1, c0, c1);

  fmpz_t d0;
  fmpz_init(d0);
  fmpz_add(d0, dabit00->bp, a0);  // x0 + a0
  fmpz_t e0;
  fmpz_init(e0);
  fmpz_add(e0, dabit10->bp, b0);  // y0 + b0
  fmpz_t d1;
  fmpz_init(d1);
  fmpz_add(d1, dabit01->bp, a1);  // x1 + a1
  fmpz_t e1;
  fmpz_init(e1);
  fmpz_add(e1, dabit11->bp, b1);  // y1 + b1

  fmpz_t d, e;
  fmpz_init(d);
  fmpz_init(e);
  fmpz_add(d, d0, d1);  // x - a
  fmpz_mod(d, d, Int_Modulus);
  fmpz_add(e, e0, e1);  // y - b
  fmpz_mod(e, e, Int_Modulus);

  // [xy] = [c] + [x] e + [y] d - de
  fmpz_set(prod0, c0);
  fmpz_addmul(prod0, dabit00->bp, e);
  fmpz_addmul(prod0, dabit10->bp, d);
  fmpz_submul(prod0, d, e);
  fmpz_mod(prod0, prod0, Int_Modulus);

  fmpz_set(prod1, c1);
  fmpz_addmul(prod1, dabit01->bp, e);
  fmpz_addmul(prod1, dabit11->bp, d);
  fmpz_mod(prod1, prod1, Int_Modulus);

  fmpz_add(dabit0->bp, dabit00->bp, dabit10->bp);
  fmpz_submul_ui(dabit0->bp, prod0, 2);
  fmpz_mod(dabit0->bp, dabit0->bp, Int_Modulus);

  fmpz_add(dabit1->bp, dabit01->bp, dabit11->bp);
  fmpz_submul_ui(dabit1->bp, prod1, 2);
  fmpz_mod(dabit1->bp, dabit1->bp, Int_Modulus);
}

void localTest(const size_t n) {
  std::cout << "Running Local Test" << std::endl;

  fmpz_t tmp;
  fmpz_init(tmp);

  EdaBit* ebit00 = new EdaBit(n);  // 0's share of 0
  EdaBit* ebit01 = new EdaBit(n);  // 1's share of 0
  makeLocalEdaBit(ebit00, ebit01, n);

  EdaBit* ebit10 = new EdaBit(n);  // 0's share of 1
  EdaBit* ebit11 = new EdaBit(n);  // 1's share of 1
  makeLocalEdaBit(ebit10, ebit11, n);

  // 0 has ebit00, ebit01. 1 has ebit10, ebit11
  // They move things around
  // 0 has ebit00, ebit10. 1 has ebit01, ebit11

  // Final bits
  EdaBit* ebit0 = new EdaBit(n);
  EdaBit* ebit1 = new EdaBit(n);

  // r'
  fmpz_add(ebit0->r, ebit00->r, ebit10->r);
  fmpz_mod(ebit0->r, ebit0->r, Int_Modulus);

  fmpz_add(ebit1->r, ebit01->r, ebit11->r);
  fmpz_mod(ebit1->r, ebit1->r, Int_Modulus);

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
    a0 = fmpz_tstbit(tmp, 0);
    a1 = fmpz_tstbit(tmp, 1);
    b0 = fmpz_tstbit(tmp, 2);
    b1 = fmpz_tstbit(tmp, 3);
    c0 = fmpz_tstbit(tmp, 4);
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
  DaBit* dabit0 = new DaBit();
  DaBit* dabit1 = new DaBit();
  
  makeSharedDabitLocally(dabit0, dabit1);

  std::cout << "conversion dabit0\n"; dabit0->print();
  std::cout << "conversion dabit1\n"; dabit1->print();

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
  fmpz_mul_ui(tmp, bp0, 1ul << n);
  fmpz_sub(ebit0->r, ebit0->r, tmp);
  fmpz_mod(ebit0->r, ebit0->r, Int_Modulus);

  fmpz_mul_ui(tmp, bp1, 1ul << n);
  fmpz_sub(ebit1->r, ebit1->r, tmp);
  fmpz_mod(ebit1->r, ebit1->r, Int_Modulus);

  std::cout << "Final ebit0:\n"; ebit0->print();
  std::cout << "Final ebit1:\n"; ebit1->print();

  fmpz_add(tmp, ebit0->r, ebit1->r);
  fmpz_mod(tmp, tmp, Int_Modulus);
  std::cout << "Validation" << std::endl;
  std::cout << " sum = "; fmpz_print(tmp); std::cout << std::endl;
  fmpz_t tmp2;
  fmpz_init(tmp2);
  fmpz_from_bool_array(tmp, ebit0->b, n);
  fmpz_from_bool_array(tmp2, ebit1->b, n);
  fmpz_xor(tmp, tmp, tmp2);
  fmpz_clear(tmp2);
  std::cout << " xor = "; fmpz_print(tmp); std::cout << std::endl;
}

void runServerTest(const int server_num, const int serverfd, const size_t n) {

  BooleanBeaverTriple* triples = gen_boolean_beaver_triples(server_num, n);

  // TODO: use proper triple generation
  BeaverTriple* triple = new BeaverTriple();
  if (server_num == 0) {
    BeaverTriple* other_triple = new BeaverTriple();
    makeBeaverTripleLocally(triple->A, other_triple->A,
                            triple->B, other_triple->B,
                            triple->C, other_triple->C);
    send_BeaverTriple(serverfd, other_triple);
    delete other_triple;
  } else {
    recv_BeaverTriple(serverfd, triple);
  }

  DaBit* dabit = generateDaBit(serverfd, server_num, triple);

  delete triple;

  EdaBit* ebit = generateEdaBit(serverfd, server_num, n, triples, dabit);

  delete[] triples;
  delete dabit;

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
    fmpz_t tmp2;
    fmpz_init(tmp2);
    fmpz_from_bool_array(tmp, ebit->b, n);
    fmpz_from_bool_array(tmp2, ebit_other->b, n);
    fmpz_xor(tmp, tmp, tmp2);
    fmpz_clear(tmp2);
    std::cout << " xor = "; fmpz_print(tmp); std::cout << std::endl;

    fmpz_clear(tmp);
    delete ebit_other;
  }

  sleep(1);

  fmpz_t x;
  fmpz_init(x);
  fmpz_t xp;
  fmpz_init(xp);

  // Use bit
  std::cout << "converting with bit" << std::endl;
  if (server_num == 0) {
    fmpz_t pow;
    fmpz_init_set_si(pow, 1l << n);
    fmpz_t x_true;  // true value
    fmpz_init(x_true);
    fmpz_randm(x_true, seed, pow);
    fmpz_clear(pow);

    fmpz_randm(x, seed, pow);
    fmpz_t x_other;
    fmpz_init(x_other);
    fmpz_xor(x_other, x_true, x);
    std::cout << "True: "; fmpz_print(x_true);
    std::cout << " = "; fmpz_print(x);
    std::cout << " ^ "; fmpz_print(x_other);
    std::cout << std::endl;
    send_fmpz(serverfd, x_other);
    fmpz_clear(x_other);

    fmpz_randm(xp, seed, Int_Modulus);
    fmpz_t xp_other;
    fmpz_init(xp_other);
    fmpz_sub(xp_other, x_true, xp);
    fmpz_mod(xp_other, xp_other, Int_Modulus);
    std::cout << " = "; fmpz_print(xp);
    std::cout << " + "; fmpz_print(xp_other);
    std::cout << std::endl;
    send_fmpz(serverfd, xp_other);
    fmpz_clear(xp_other);
  } else {
    recv_fmpz(serverfd, x);
    recv_fmpz(serverfd, xp);
  }

  BooleanBeaverTriple* newtriples = gen_boolean_beaver_triples(server_num, n);

  bool valid = validate_shares_match(serverfd, server_num, x, xp, n, ebit, newtriples);

  std::cout << "Server " << server_num << " shares match: " << std::boolalpha << valid << std::endl;
}

void serverTest(const size_t n) {
  std::cout << "Running server test" << std::endl;
  int sockfd = init_receiver();

  pid_t pid = fork();
  if (pid == 0) {
    int cli_sockfd = init_sender();

    fmpz_t tmp;
    fmpz_init(tmp);
    fmpz_clear(tmp);

    runServerTest(0, cli_sockfd, n);
    close(cli_sockfd);
  } else if (pid > 0) {
    int newsockfd = accept_receiver(sockfd);

    // alter randomness to be different from the sender
    fmpz_t tmp;
    fmpz_init(tmp);
    size_t rand_adjust = 4;
    for (int i = 0; i < rand_adjust; i++)
      fmpz_randm(tmp, seed, Int_Modulus);
    fmpz_clear(tmp);

    runServerTest(1, newsockfd, n);

    close(newsockfd);
    close(sockfd);
  }
}

int main(int argc, char* argv[]) {
  init_constants();
  
  // Should be < log_2 Int_Modulus / 2
  size_t n = 32;
  std::cout << "n = " << n << std::endl;

  // random adjusting. different numbers adjust seed.
  fmpz_t tmp;
  fmpz_init(tmp);
  size_t rand_adjust = 5;
  for (int i = 0; i < rand_adjust; i++)
    fmpz_randm(tmp, seed, Int_Modulus);
  fmpz_clear(tmp);

  localTest(n);
  std::cout << std::endl;
  serverTest(n);

  clear_constants();
  return 0;
}