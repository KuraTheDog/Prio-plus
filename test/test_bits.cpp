#undef NDEBUG
#include <assert.h>
#include <iostream>

#include "utils_test_connect.h"
#include "../constants.h"
#include "../correlated.h"
#include "../fmpz_utils.h"
#include "../net_share.h"

const size_t batch_size = 100; // flexible
const size_t N = 80;           // Must be >= 2

const size_t num_bits = 3;     // Must be >= 3
const bool lazy = false;

void test_multiplyBoolShares(const size_t N, const int server_num, const int serverfd,
    PrecomputeStore* store) {
  bool* const x = new bool[N];
  bool* const y = new bool[N];
  bool* actual = new bool[N];
  memset(x, 0, N * sizeof(bool));
  memset(y, 0, N * sizeof(bool));
  for (unsigned int i = 0; i < N; i++) {
    // Random base
    // x[i] = (i>>3)%2;
    // y[i] = (i>>4)%2;
    // Actuals
    if (server_num == 0) {
      x[i] ^= i%2;
      y[i] ^= (i>>1)%2;
      actual[i] = ((i>>1)%2) && (i%2);
    }
  }

  bool* z = new bool[N];
  store->multiply_BoolShares(N, x, y, z);

  if (server_num == 0) {
    bool* z_other = new bool[N];
    recv_bool_batch(serverfd, z_other, N);
    for (unsigned int i = 0; i < N; i++) {
      assert((z[i] ^ z_other[i]) == actual[i]);
    }
    delete[] z_other;
  } else {
    send_bool_batch(serverfd, z, N);
  }
  delete[] x;
  delete[] y;
  delete[] z;
  delete[] actual;
}

void test_multiplyBoolShares_cross(const int server_num, const int serverfd,
    PrecomputeStore* store) {
  const size_t N = 2;
  const size_t a = 2;
  const size_t b = 3;
  bool* const x = new bool[N * a];
  bool* const y = new bool[N * b];
  memset(x, 0, N * a * sizeof(bool));
  memset(y, 0, N * b * sizeof(bool));
  // (1, 0) x (1, 1, 0)
  // (0, 1) x (1, 0, 1)
  if (server_num == 0) {
    x[0] = true; x[1] = false; x[2] = true; x[3] = false;
    y[0] = false; y[1] = true; y[2] = false; y[3] = true; y[4] = true; y[5] = false;
  } else {
    x[0] = false; x[1] = false; x[2] = true; x[3] = true;
    y[0] = true; y[1] = false; y[2] = false; y[3] = false; y[4] = true; y[5] = true;
  }
  bool* z = new bool[N * a * b];

  store->multiply_BoolShares_cross(N, a, b, x, y, z);

  reveal_bool_batch(serverfd, x, N * a);
  reveal_bool_batch(serverfd, y, N * b);
  reveal_bool_batch(serverfd, z, N * a * b);
  if (server_num == 0) {
    for (unsigned int i = 0; i < N; i++) {
      for (unsigned int j = 0; j < a; j++) {
        for (unsigned int k = 0; k < b; k++) {
          int idx = i * a * b + (j * b + k);
          int x_idx = i * a + j;
          int y_idx = i * b + k;
          // std::cout << "z[" << idx << "] (" << i << ", " << j << ", " << k << ") = " << z[idx];
          // std::cout << " = " << x[x_idx] << " * " << y[y_idx] << std::endl;
          assert(z[idx] == x[x_idx] * y[y_idx]);
        }
      }
    }
  }

  delete[] x;
  delete[] y;
  delete[] z;
}

void test_addBinaryShares(const size_t N, const size_t* const nbits,
    const int server_num, const int serverfd, PrecomputeStore* store) {
  bool** const x = new bool*[N];
  bool** const y = new bool*[N];
  bool** const z = new bool*[N];
  for (unsigned int i = 0; i < N; i++) {
    x[i] = new bool[nbits[i]];
    y[i] = new bool[nbits[i]];
    z[i] = new bool[nbits[i]];
    memset(x[i], 0, nbits[i] * sizeof(bool));
    memset(y[i], 0, nbits[i] * sizeof(bool));
    memset(z[i], 0, nbits[i] * sizeof(bool));
  }

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
  store->add_BinaryShares(N, nbits, x, y, z, carry);

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

void test_b2a_single(const size_t N, const int server_num, const int serverfd,
    ShareConverter* store) {
  bool* x = new bool[N];
  for (unsigned int i = 0; i < N; i++) {
    x[i] = (i+1)%2;  // common mask
    if (server_num == 0)
      x[i] ^= (i>>2) % 2;    // true
  }

  fmpz_t* xp; new_fmpz_array(&xp, N);
  store->b2a_single(N, x, xp);

  if (server_num == 0) {
    reveal_fmpz_batch(serverfd, xp, N);

    for (unsigned int i = 0; i < N; i++)
      assert(fmpz_equal_ui(xp[i], (i>>2) % 2));
  } else {
    reveal_fmpz_batch(serverfd, xp, N);
  }

  delete[] x;

  clear_fmpz_array(xp, N);
}

void test_b2a_multi(const size_t N, const int server_num, const int serverfd,
    ShareConverter* store) {
  uint64_t x[N];
  x[0] = server_num == 0 ? 3 : 6;
  x[1] = server_num == 0 ? 5 : 4;

  fmpz_t* xp; new_fmpz_array(&xp, N);
  store->b2a_multi(N, num_bits, x, xp);

  if (server_num == 0) {
    fmpz_t tmp; fmpz_init(tmp);

    recv_fmpz(serverfd, tmp);
    fmpz_mod_add(tmp, tmp, xp[0], mod_ctx);
    assert(fmpz_equal_ui(tmp, 5));  // 3 ^ 6 = 5

    recv_fmpz(serverfd, tmp);
    fmpz_mod_add(tmp, tmp, xp[1], mod_ctx);
    assert(fmpz_equal_ui(tmp, 1));  // 5 ^ 4 = 1

    fmpz_clear(tmp);
  } else {
    send_fmpz(serverfd, xp[0]);
    send_fmpz(serverfd, xp[1]);
  }

  clear_fmpz_array(xp, N);
}

void runServerTest(const int server_num, const int serverfd) {
  OT_Wrapper* ot0 = new OT_Wrapper(server_num == 0 ? nullptr : "127.0.0.1", 60051);
  OT_Wrapper* ot1 = new OT_Wrapper(server_num == 1 ? nullptr : "127.0.0.1", 60052);
  PrecomputeStore* store_pre = new PrecomputeStore(serverfd, server_num, ot0, ot1, batch_size, lazy);
  OTCorrelatedStore* store_ot = new OTCorrelatedStore(serverfd, server_num, ot0, ot1);
  store_pre->maybe_Update();

  size_t* bits_arr = new size_t[N];
  for (unsigned int i = 0; i < N; i++)
    bits_arr[i] = (i % 2 == 0 ? num_bits : 2 * num_bits);
  bits_arr[1] = num_bits * 2;

  auto start = clock_start();

  for (int i = 0; i < 1; i++) {
    std::cout << "iteration: " << i << std::endl;
    start = clock_start();

    test_multiplyBoolShares(N, server_num, serverfd, store_pre);
    std::cout << "mul bool timing : " << sec_from(start) << std::endl; start = clock_start();

    test_multiplyBoolShares_cross(server_num, serverfd, store_pre);
    std::cout << "mul bool cross timing : " << sec_from(start) << std::endl; start = clock_start();

    test_addBinaryShares(N, bits_arr, server_num, serverfd, store_pre);
    std::cout << "add bin timing : " << sec_from(start) << std::endl; start = clock_start();

    test_b2a_single(N, server_num, serverfd, store_pre);
    std::cout << "b2a da single timing : " << sec_from(start) << std::endl; start = clock_start();

    test_b2a_multi(N, server_num, serverfd, store_pre);
    std::cout << "b2a da multi timing : " << sec_from(start) << std::endl; start = clock_start();

    test_b2a_single(N, server_num, serverfd, store_ot);
    std::cout << "b2a ot single timing : " << sec_from(start) << std::endl; start = clock_start();

    test_b2a_multi(N, server_num, serverfd, store_ot);
    std::cout << "b2a ot multi timing : " << sec_from(start) << std::endl; start = clock_start();
  }

  store_pre->print_Sizes();

  delete ot0;
  delete ot1;
  delete[] bits_arr;
  delete store_pre;
  delete store_ot;
}

void serverTest() {
  std::cout << "Running server test" << std::endl;

  std::thread t0([&]() {
    int cli_sockfd = init_sender();

    runServerTest(0, cli_sockfd);

    close(cli_sockfd);
  });

  std::thread t1([&]() {
    int sockfd = init_receiver();
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
    close(sockfd);
  });

  t0.join();
  t1.join();
}

int main(int argc, char* argv[]) {
  init_constants();

  int server_num = -1;
  if (argc >= 2){
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
