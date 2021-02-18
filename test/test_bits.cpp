#include <iostream>

#include "utils_test_connect.h"
#include "../constants.h"
#include "../edabit.h"
#include "../fmpz_utils.h"
#include "../net_share.h"

void test_multiplyBoolShares(const size_t N, const int server_num, const int serverfd, CorrelatedStore* store) {
  bool* x = new bool[N];
  bool* y = new bool[N];
  if (server_num == 0) {
    x[0] = true; x[1] = false;
    y[0] = false; y[1] = true;
  } else {
    x[0] = true; x[1] = true;
    x[1] = true; y[1] = false;
  }

  bool* z = store->multiplyBoolShares(N, x, y);

  if (server_num == 0) {
    bool z_other;
    recv_bool(serverfd, z_other);
    assert((z[0] ^ z_other) == false);  // (1 ^ 1) (0 ^ 1) = 0
    recv_bool(serverfd, z_other);
    assert((z[1] ^ z_other) == true);   // (0 ^ 1) (1 ^ 0) = 1
  } else {
    send_bool(serverfd, z[0]);
    send_bool(serverfd, z[1]);
  }
  delete[] x;
  delete[] y;
  delete[] z;
}

void test_multiplyArithmeticShares(const size_t N, const int server_num, const int serverfd, CorrelatedStore* store) {
  fmpz_t* x; new_fmpz_array(&x, N);
  fmpz_t* y; new_fmpz_array(&y, N);

  if (server_num == 0) {
    fmpz_set_ui(x[0], 37);
    fmpz_set_ui(x[1], 2);
    fmpz_set_si(y[0], -84); fmpz_mod(y[0], y[0], Int_Modulus);
    fmpz_set_ui(y[1], 69);
  } else {
    fmpz_set_si(x[0], -30); fmpz_mod(x[0], x[0], Int_Modulus);
    fmpz_set_ui(x[1], 11);
    fmpz_set_ui(y[0], 90);
    fmpz_set_si(y[1], -62); fmpz_mod(y[1], y[1], Int_Modulus);
  }

  fmpz_t* z = store->multiplyArithmeticShares(N, x, y);

  if (server_num == 0) {
    fmpz_t tmp; fmpz_init(tmp);

    recv_fmpz(serverfd, tmp);
    fmpz_add(tmp, tmp, z[0]);
    fmpz_mod(tmp, tmp, Int_Modulus);
    assert(fmpz_equal_ui(tmp, 42)); // 7 * 6 = 42

    recv_fmpz(serverfd, tmp);
    fmpz_add(tmp, tmp, z[1]);
    fmpz_mod(tmp, tmp, Int_Modulus);
    assert(fmpz_equal_ui(tmp, 91)); // 13 * 7 = 91

    fmpz_clear(tmp);
  } else {
    send_fmpz(serverfd, z[0]);
    send_fmpz(serverfd, z[1]);
  }

  clear_fmpz_array(x, N);
  clear_fmpz_array(y, N);
  clear_fmpz_array(z, N);
}

void test_addBinaryShares(const size_t N, const int server_num, const int serverfd, CorrelatedStore* store) {
  bool** x = new bool*[N];
  bool** y = new bool*[N];
  bool** z = new bool*[N];
  for (unsigned int i = 0; i < N; i++) {
    x[i] = new bool[3];
    y[i] = new bool[3];
    z[i] = new bool[3];
  }

  if (server_num == 0) {
    x[0][0] = 1; x[0][1] = 0; x[0][2] = 1;  // 3
    y[0][0] = 1; y[0][1] = 1; y[0][2] = 0;  // 2

    x[1][0] = 0; x[1][1] = 0; x[1][2] = 1;  // 5
    y[1][0] = 1; y[1][1] = 0; y[1][2] = 0;  // 4
  } else {
    x[0][0] = 0; x[0][1] = 1; x[0][2] = 1;
    y[0][0] = 1; y[0][1] = 0; y[0][2] = 0;

    x[1][0] = 1; x[1][1] = 0; x[1][2] = 0;
    y[1][0] = 1; y[1][1] = 0; y[1][2] = 1;
  }

  bool* carry = store->addBinaryShares(N, x, y, z);

  if (server_num == 0) {
    bool other;
    // 3 + 2 = 5
    recv_bool(serverfd, other);
    assert((z[0][0] ^ other) == true);
    recv_bool(serverfd, other);
    assert((z[0][1] ^ other) == false);
    recv_bool(serverfd, other);
    assert((z[0][2] ^ other) == true);
    recv_bool(serverfd, other);
    assert((carry[0] ^ other) == false);
    // 4 + 5 = 9
    recv_bool(serverfd, other);
    assert((z[1][0] ^ other) == true);
    recv_bool(serverfd, other);
    assert((z[1][1] ^ other) == false);
    recv_bool(serverfd, other);
    assert((z[1][2] ^ other) == false);
    recv_bool(serverfd, other);
    assert((carry[1] ^ other) == true);
  } else {
    send_bool(serverfd, z[0][0]);
    send_bool(serverfd, z[0][1]);
    send_bool(serverfd, z[0][2]);
    send_bool(serverfd, carry[0]);

    send_bool(serverfd, z[1][0]);
    send_bool(serverfd, z[1][1]);
    send_bool(serverfd, z[1][2]);
    send_bool(serverfd, carry[1]);

  }

  for (unsigned int i = 0; i < N; i++) {
    delete[] x[i];
    delete[] y[i];
    delete[] z[i];
  }
  delete[] x;
  delete[] y;
  delete[] z;
  delete[] carry;
}

void test_b2a_daBit(const size_t N, const int server_num, const int serverfd, CorrelatedStore* store) {
  bool* x = new bool[N];
  if (server_num == 0) {
    x[0] = true; x[1] = false;
  } else {
    x[0] = true; x[1] = true;
  }

  fmpz_t* xp = store->b2a_daBit(N, x);

  if (server_num == 0) {
    fmpz_t tmp; fmpz_init(tmp);

    recv_fmpz(serverfd, tmp);
    fmpz_add(tmp, tmp, xp[0]);
    fmpz_mod(tmp, tmp, Int_Modulus);
    assert(fmpz_equal_ui(tmp, 0));  // 1 ^ 1 = 0

    recv_fmpz(serverfd, tmp);
    fmpz_add(tmp, tmp, xp[1]);
    fmpz_mod(tmp, tmp, Int_Modulus);
    assert(fmpz_equal_ui(tmp, 1));  // 0 ^ 1 = 1

    fmpz_clear(tmp);
  } else {
    send_fmpz(serverfd, xp[0]);
    send_fmpz(serverfd, xp[1]);
  }

  delete[] x;

  clear_fmpz_array(xp, N);
} 

void test_b2a_edaBit(const size_t N, const int server_num, const int serverfd, CorrelatedStore* store) {
  fmpz_t* x; new_fmpz_array(&x, N);

  if (server_num == 0) {
    fmpz_set_ui(x[0], 3);
    fmpz_set_ui(x[1], 5);
  } else {
    fmpz_set_ui(x[0], 6);
    fmpz_set_ui(x[1], 4);
  }

  fmpz_t* xp = store->b2a_edaBit(N, x);

  if (server_num == 0) {
    fmpz_t tmp; fmpz_init(tmp);

    recv_fmpz(serverfd, tmp);
    fmpz_add(tmp, tmp, xp[0]);
    fmpz_mod(tmp, tmp, Int_Modulus);
    assert(fmpz_equal_ui(tmp, 5));  // 3 ^ 6 = 5

    recv_fmpz(serverfd, tmp);
    fmpz_add(tmp, tmp, xp[1]);
    fmpz_mod(tmp, tmp, Int_Modulus);
    assert(fmpz_equal_ui(tmp, 1));  // 5 ^ 4 = 1

    fmpz_clear(tmp);
  } else {
    send_fmpz(serverfd, xp[0]);
    send_fmpz(serverfd, xp[1]);
  }

  clear_fmpz_array(x, N);
  clear_fmpz_array(xp, N);
}

void test_validateSharesMatch(const size_t N, const int server_num, const int serverfd, CorrelatedStore* store) {
  fmpz_t* x2; new_fmpz_array(&x2, N);
  fmpz_t* xp; new_fmpz_array(&xp, N);

  if (server_num == 0) {
    fmpz_set_ui(x2[0], 3);
    fmpz_set_ui(x2[1], 5);
    fmpz_set_ui(xp[0], 42);
    fmpz_set_si(xp[1], -68); fmpz_mod(xp[1], xp[1], Int_Modulus);
  } else {
    fmpz_set_ui(x2[0], 6);
    fmpz_set_ui(x2[1], 4);
    fmpz_set_si(xp[0], -37); fmpz_mod(xp[0], xp[0], Int_Modulus);
    fmpz_set_ui(xp[1], 69);
  }

  bool* ret = store->validateSharesMatch(N, x2, xp);

  assert(ret[0]);
  assert(ret[1]);

  delete[] ret;
  clear_fmpz_array(x2, N);
  clear_fmpz_array(xp, N);
}

void runServerTest(const int server_num, const int serverfd) {
  const bool lazy = false;
  const size_t nbits = 3;
  const size_t batch_size = 50000;
  const size_t N = 50000;  // >= 2
  CorrelatedStore* store = new CorrelatedStore(serverfd, server_num, "127.0.0.1", "127.0.0.1", nbits, batch_size, lazy);

  store->maybeUpdate();

  std::cout << std::endl;

  for (int i = 0; i < 1; i++) {
    std::cout << "iteration: " << i << std::endl;
    std::cout << "mul bool" << std::endl;
    test_multiplyBoolShares(N, server_num, serverfd, store);
    std::cout << "mul arith" << std::endl;
    test_multiplyArithmeticShares(N, server_num, serverfd, store);
    std::cout << "add bin" << std::endl;
    test_addBinaryShares(N, server_num, serverfd, store);
    std::cout << "b2a da" << std::endl;
    test_b2a_daBit(N, server_num, serverfd, store);
    std::cout << "b2a ed" << std::endl;
    test_b2a_edaBit(N, server_num, serverfd, store);
    std::cout << "validate" << std::endl;
    test_validateSharesMatch(N, server_num, serverfd, store);
  }

  delete store;
}

void serverTest() {
  std::cout << "Running server test" << std::endl;
  int sockfd = init_receiver();

  pid_t pid = fork();
  if (pid == 0) {
    int cli_sockfd = init_sender();

    runServerTest(0, cli_sockfd);
    close(cli_sockfd);
  } else {
    int newsockfd = accept_receiver(sockfd);

    // alter randomness to be different from the sender
    fmpz_t tmp;
    fmpz_init(tmp);
    size_t rand_adjust = 4;
    for (unsigned int i = 0; i < rand_adjust; i++)
      fmpz_randm(tmp, seed, Int_Modulus);
    fmpz_clear(tmp);

    runServerTest(1, newsockfd);

    close(newsockfd);
  }

  close(sockfd);
}

int main(int argc, char* argv[]) {
  init_constants();

  int server_num = -1;
  if(argc >= 2){
    server_num = atoi(argv[1]);
  }

  // random adjusting. different numbers adjust seed.
  fmpz_t tmp;
  fmpz_init(tmp);
  size_t rand_adjust = 5;
  for (unsigned int i = 0; i < rand_adjust; i++)
    fmpz_randm(tmp, seed, Int_Modulus);
  fmpz_clear(tmp);

  if (server_num == -1) {
    serverTest();
  } else if (server_num == 0) {
    int sockfd = init_receiver();
    int newsockfd = accept_receiver(sockfd);
    runServerTest(0, newsockfd);
    close(newsockfd);
    close(sockfd);
  } else if (server_num == 1) {
    int cli_sockfd = init_sender();
    runServerTest(1, cli_sockfd);
    close(cli_sockfd);
  }

  std::cout << "All tests passed" << std::endl;

  clear_constants();
  return 0;
}