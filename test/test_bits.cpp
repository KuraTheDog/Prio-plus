#undef NDEBUG
#include <assert.h>
#include <iostream>

#include "utils_test_connect.h"
#include "../constants.h"
#include "../correlated.h"
#include "../fmpz_utils.h"
#include "../net_share.h"

const size_t batch_size = 100; // flexible
const size_t N = 8;           // Must be >= 2

const size_t num_bits = 3;     // Must be >= 3
const bool do_fork = true;     // false to do leaks testing
const bool lazy = false;

void test_multiplyBoolShares(const size_t N, const int server_num,
    const int serverfd, CorrelatedStore* store) {
  bool* const x = new bool[N];
  bool* const y = new bool[N];
  if (server_num == 0) {
    x[0] = true; x[1] = false;
    y[0] = false; y[1] = true;
  } else {
    x[0] = true; x[1] = true;
    x[1] = true; y[1] = false;
  }

  bool* z = new bool[N];
  store->multiplyBoolShares(N, x, y, z);

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

void test_multiplyArithmeticShares(const size_t N, const int server_num,
    const int serverfd, CorrelatedStore* store) {
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

  fmpz_t* z; new_fmpz_array(&z, N);
  store->multiplyArithmeticShares(N, x, y, z);

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

void test_multiplyBoolArith(const size_t N, const int server_num,
    const int serverfd, CorrelatedStore* store) {
  bool* const b = new bool[4];
  fmpz_t* x; new_fmpz_array(&x, N * 4);

  // bit value: 0, 1, 1, 0
  if (server_num == 0) {
    b[0] = false; b[1] = false; b[2] = true; b[3] = true;
    for (unsigned int i = 0; i < N*4; i++) {
      fmpz_set_si(x[i], i+1);
    }
  } else {
    b[0] = false; b[1] = true; b[2] = false; b[3] = true;
    for (unsigned int i = 0; i < N*4; i++) {
      fmpz_set_si(x[i], (i+1)*1000);
    }
  }

  fmpz_t* z; new_fmpz_array(&z, N * 4);
  store->multiplyBoolArithFlat(N, 4, b, x, z);

  if (server_num == 0) {
    fmpz_t* z_other; new_fmpz_array(&z_other, N * 4);
    recv_fmpz_batch(serverfd, z_other, N * 4);
    for (unsigned int i = 0; i < N * 4; i++) {
      fmpz_add(z[i], z[i], z_other[i]);
      fmpz_mod(z[i], z[i], Int_Modulus);

      // std::cout << "z[" << i << "] = "; fmpz_print(z[i]); std::cout << std::endl;
    }
    clear_fmpz_array(z_other, N);

    for (unsigned int j = 0; j < 4; j++) {
      // std::cout << "bit: " << j << std::endl;
      for (unsigned int i = 0; i < N; i++) {
        const size_t idx = i * 4 + j;
        // std::cout << "N = " << i << ", idx = " << idx << ": "; fmpz_print(z[idx]); std::cout << std::endl;
        if (j == 0 or j == 3) {
          assert(fmpz_equal_ui(z[idx], 0));
        } else {
          assert(fmpz_equal_ui(z[idx], (idx+1) * 1001));
        }
      }
    }
  } else {
    send_fmpz_batch(serverfd, z, N * 4);
  }

  delete[] b;
  clear_fmpz_array(x, N * 4);
  clear_fmpz_array(z, N * 4);
}

void test_addBinaryShares(const size_t N, const size_t* const nbits, const int server_num,
    const int serverfd, CorrelatedStore* store) {
  bool** const x = new bool*[N];
  bool** const y = new bool*[N];
  bool** const z = new bool*[N];
  for (unsigned int i = 0; i < N; i++) {
    x[i] = new bool[nbits[i]];
    y[i] = new bool[nbits[i]];
    z[i] = new bool[nbits[i]];
  }
  memset(x[0], 0, nbits[0]);
  memset(x[1], 0, nbits[1]);
  memset(y[0], 0, nbits[0]);
  memset(y[1], 0, nbits[1]);

  if (server_num == 0) {
    x[0][0] = 1; x[0][1] = 0; x[0][2] = 1;  // 3
    y[0][0] = 1; y[0][1] = 1; y[0][2] = 0;  // 2

    x[1][0] = 0; x[1][1] = 0; x[1][2] = 1;  x[1][3] = 0; // 13
    y[1][0] = 1; y[1][1] = 0; y[1][2] = 0;  y[1][3] = 1; // 12
  } else {
    x[0][0] = 0; x[0][1] = 1; x[0][2] = 1;
    y[0][0] = 1; y[0][1] = 0; y[0][2] = 0;

    x[1][0] = 1; x[1][1] = 0; x[1][2] = 0; x[1][3] = 1;
    y[1][0] = 1; y[1][1] = 0; y[1][2] = 1; y[1][3] = 0;
  }

  bool* carry = new bool[N];
  store->addBinaryShares(N, nbits, x, y, z, carry);

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
    // 13 + 12 = 25
    recv_bool(serverfd, other);
    assert((z[1][0] ^ other) == true);
    recv_bool(serverfd, other);
    assert((z[1][1] ^ other) == false);
    recv_bool(serverfd, other);
    assert((z[1][2] ^ other) == false);
    recv_bool(serverfd, other);
    assert((z[1][3] ^ other) == true);
    recv_bool(serverfd, other);
    assert((z[1][4] ^ other) == true);
    recv_bool(serverfd, other);
    assert((z[1][5] ^ other) == false);
    recv_bool(serverfd, other);
    assert((carry[1] ^ other) == false);
  } else {
    send_bool(serverfd, z[0][0]);
    send_bool(serverfd, z[0][1]);
    send_bool(serverfd, z[0][2]);
    send_bool(serverfd, carry[0]);

    send_bool(serverfd, z[1][0]);
    send_bool(serverfd, z[1][1]);
    send_bool(serverfd, z[1][2]);
    send_bool(serverfd, z[1][3]);
    send_bool(serverfd, z[1][4]);
    send_bool(serverfd, z[1][5]);
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

void test_b2a_daBit_single(const size_t N, const int server_num,
    const int serverfd, CorrelatedStore* store) {
  bool* x = new bool[N];
  if (server_num == 0) {
    x[0] = true; x[1] = false;
  } else {
    x[0] = true; x[1] = true;
  }

  fmpz_t* xp; new_fmpz_array(&xp, N);
  store->b2a_daBit_single(N, x, xp);

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

void test_b2a_multi(const size_t N, const size_t* const nbits, const int server_num,
    const int serverfd, CorrelatedStore* store) {
  fmpz_t* x; new_fmpz_array(&x, N);

  if (server_num == 0) {
    fmpz_set_ui(x[0], 3);
    fmpz_set_ui(x[1], 5);
  } else {
    fmpz_set_ui(x[0], 6);
    fmpz_set_ui(x[1], 4);
  }

  fmpz_t* xp; new_fmpz_array(&xp, N);
  store->b2a_daBit_multi(N, nbits, x, xp);

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

void test_b2a_ot(const size_t N, const size_t* const nbits, const int server_num,
    const int serverfd, CorrelatedStore* store) {
  fmpz_t* x; new_fmpz_array(&x, N);

  if (server_num == 0) {
    fmpz_set_ui(x[0], 3);
    fmpz_set_ui(x[1], 5);
  } else {
    fmpz_set_ui(x[0], 6);
    fmpz_set_ui(x[1], 4);
  }

  const size_t mod = fmpz_get_ui(Int_Modulus);
  fmpz_t* xp; new_fmpz_array(&xp, N);
  store->b2a_ot(N, 1, nbits, x, xp, mod);

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

void runServerTest(const int server_num, const int serverfd) {
  OT_Wrapper* ot0 = new OT_Wrapper(server_num == 0 ? nullptr : "127.0.0.1", 60051);
  OT_Wrapper* ot1 = new OT_Wrapper(server_num == 1 ? nullptr : "127.0.0.1", 60052);
  CorrelatedStore* store = new CorrelatedStore(serverfd, server_num, ot0, ot1, batch_size, lazy, do_fork);

  store->maybeUpdate();

  std::cout << std::endl;

  size_t* bits_arr = new size_t[N];
  for (unsigned int i = 0; i < N; i++)
    bits_arr[i] = (i % 2 == 0 ? num_bits : 2 * num_bits);
  bits_arr[1] = num_bits * 2;

  auto start = clock_start();

  for (int i = 0; i < 1; i++) {
    std::cout << "iteration: " << i << std::endl;
    start = clock_start();

    /* Unused
    test_multiplyBoolShares(N, server_num, serverfd, store);
    std::cout << "mul bool timing : " << sec_from(start) << std::endl; start = clock_start();
    */

    test_multiplyArithmeticShares(N, server_num, serverfd, store);
    std::cout << "mul arith timing : " << sec_from(start) << std::endl; start = clock_start();

    test_multiplyBoolArith(N, server_num, serverfd, store);
    std::cout << "mul bool arith timing : " << sec_from(start) << std::endl; start = clock_start();

    /* Unused
    test_addBinaryShares(N, bits_arr, server_num, serverfd, store);
    std::cout << "add bin timing : " << sec_from(start) << std::endl; start = clock_start();
    */

    test_b2a_daBit_single(N, server_num, serverfd, store);
    std::cout << "b2a da single timing : " << sec_from(start) << std::endl; start = clock_start();

    test_b2a_multi(N, bits_arr, server_num, serverfd, store);
    std::cout << "b2a da multi timing : " << sec_from(start) << std::endl; start = clock_start();

    test_b2a_ot(N, bits_arr, server_num, serverfd, store);
    std::cout << "b2a ot timing : " << sec_from(start) << std::endl; start = clock_start();
  }

  store->printSizes();

  delete ot0;
  delete ot1;
  delete[] bits_arr;
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
