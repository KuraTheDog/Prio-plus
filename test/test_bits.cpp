#include <iostream>

#include "utils_test_connect.h"
#include "../constants.h"
#include "../edabit.h"
#include "../fmpz_utils.h"
#include "../net_share.h"

void test_multiplyBoolShares(const int server_num, const int serverfd, CorrelatedStore* store) {
  bool* x = new bool[2];
  bool* y = new bool[2];
  if (server_num == 0) {
    x[0] = true; x[1] = false;
    y[0] = false; y[1] = true;
  } else {
    x[0] = true; x[1] = true;
    x[1] = true; y[1] = false;
  }

  bool* z = store->multiplyBoolShares(2, x, y);

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

void test_multiplyArithmeticShares(const int server_num, const int serverfd, CorrelatedStore* store) {
  fmpz_t* x; new_fmpz_array(&x, 2);
  fmpz_t* y; new_fmpz_array(&y, 2);

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

  fmpz_t* z = store->multiplyArithmeticShares(2, x, y);

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

  clear_fmpz_array(x, 2);
  clear_fmpz_array(y, 2);
  clear_fmpz_array(z, 2);
}

void test_addBinaryShares(const int server_num, const int serverfd, CorrelatedStore* store) {
  bool** x = new bool*[2];
  bool** y = new bool*[2];
  bool** z = new bool*[2];
  x[0] = new bool[3]; x[1] = new bool[3];
  y[0] = new bool[3]; y[1] = new bool[3];
  z[0] = new bool[3]; z[1] = new bool[3];

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

  bool* carry = store->addBinaryShares(2, x, y, z);

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

  delete[] x[0]; delete[] x[1];
  delete[] y[0]; delete[] y[1];
  delete[] z[0]; delete[] z[1];
  delete[] x;
  delete[] y;
  delete[] z;
  delete[] carry;
}

void test_b2a_daBit(const int server_num, const int serverfd, CorrelatedStore* store) {
  bool* x = new bool[2];
  if (server_num == 0) {
    x[0] = true; x[1] = false;
  } else {
    x[0] = true; x[1] = true;
  }

  fmpz_t* xp = store->b2a_daBit(2, x);

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

  clear_fmpz_array(xp, 2);
} 

void test_b2a_edaBit(const int server_num, const int serverfd, CorrelatedStore* store) {
  fmpz_t* x; new_fmpz_array(&x, 2);

  if (server_num == 0) {
    fmpz_set_ui(x[0], 3);
    fmpz_set_ui(x[1], 5);
  } else {
    fmpz_set_ui(x[0], 6);
    fmpz_set_ui(x[1], 4);
  }

  fmpz_t* xp = store->b2a_edaBit(2, x);

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

  clear_fmpz_array(x, 2);
  clear_fmpz_array(xp, 2);
}

void test_validateSharesMatch(const int server_num, const int serverfd, CorrelatedStore* store) {
  fmpz_t* x2; new_fmpz_array(&x2, 2);
  fmpz_t* xp; new_fmpz_array(&xp, 2);

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

  bool* ret = store->validateSharesMatch(2, x2, xp);

  assert(ret[0]);
  assert(ret[1]);

  delete[] ret;
  clear_fmpz_array(x2, 2);
  clear_fmpz_array(xp, 2);
}

void runServerTest(const int server_num, const int serverfd) {
  const bool lazy = false;
  const size_t nbits = 3;
  const size_t batch_size = 2;
  CorrelatedStore* store = new CorrelatedStore(serverfd, server_num, "127.0.0.1", "127.0.0.1", nbits, batch_size, lazy);

  for (int i = 0; i < 10; i++) {
    test_multiplyBoolShares(server_num, serverfd, store);
    test_multiplyArithmeticShares(server_num, serverfd, store);
    test_addBinaryShares(server_num, serverfd, store);
    test_b2a_daBit(server_num, serverfd, store);
    test_b2a_edaBit(server_num, serverfd, store);
    test_validateSharesMatch(server_num, serverfd, store);
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
    int cli_sockfd = init_sender();
    runServerTest(0, cli_sockfd);
    close(cli_sockfd);
  } else if (server_num == 1) {
    int sockfd = init_receiver();
    int newsockfd = accept_receiver(sockfd);
    runServerTest(1, newsockfd);
    close(newsockfd);
    close(sockfd);
  }

  std::cout << "All tests passed" << std::endl;

  clear_constants();
  return 0;
}