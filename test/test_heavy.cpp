#include <iostream>

#include "utils_test_connect.h"
#include "../constants.h"
#include "../correlated.h"
#include "../fmpz_utils.h"
#include "../net_share.h"


const size_t batch_size = 100; // flexible
const bool do_fork = true;     // false to do leaks testing
const bool lazy = false;

/*
4 bits, for each bucket/hash combo
8 inputs per, for possible randomness (0-bucket, 2 xor bits)
Tries all combinations of shares

bit = 0, all -1 into bucket 0
bit = 1, all -1 into bucket 1
bit = 2, all +1 into bucket 0
bit = 3, all +1 into bucket 1
*/

void testHeavyOT(const int server_num, const int serverfd,
                 CorrelatedStore* store) {
  const size_t nbits = 4;
  const size_t N = 8;
  const size_t n = nbits * N;

  bool* const shares_x0 = new bool[2 * n];
  bool* const shares_x1 = new bool[2 * n];
  fmpz_t* bucket0; new_fmpz_array(&bucket0, nbits);
  fmpz_t* bucket1; new_fmpz_array(&bucket1, nbits);
  bool* const valid = new bool[N]; memset(valid, 1, N);

  // Set up values for heavy convert
  for (unsigned int j = 0; j < nbits; j++) {
    const size_t bucket = j%2;
    const size_t hash = (j>>1)%2;
    // std::cout << "bit " << j << " = bucket " << bucket << ", hash " << hash << std::endl;

    for (unsigned int i = 0; i < N; i++) {
      const size_t idx = i * nbits + j;
      // std::cout << "idx(" << i << ", " << j << ") = " << idx << std::endl;
      // std::cout << " bucket" << bucket << "[" << i << "] " << (hash == 1 ? "+" : "-") << "=1" << std::endl;
      // Base
      size_t base = 0;
      if (bucket == 0) {
        base |= 1;
        base |= (hash ? 0 : 1) << 1;
        base |= ((i>>3)%2) << 3;
      } else {
        base |= ((i>>3)%2) << 1;
        base |= 1 << 2;
        base |= (hash ? 0 : 1) << 3;
      }
      // if (j % 4 == 0) {
      // std::cout << " new base: " << (base%2) << ((base>>1)%2) << ", " << ((base>>2)%2) << ((base>>3)%2) << std::endl;
      // }

      // 0: first bit same, for validity. 1 has first diff.
      size_t share = 0;
      share |= (i>>1)%2;
      share |= ((i>>2)%2) << 1;  // rand
      share |= ((i>>1)%2) << 2;  // same as first
      share |= (i%2) << 3;  // rand
      if (server_num == 1)
        share = base ^ share;

      // std::cout << "  share" << server_num << " = " << share << ": " << (share%2) << ((share>>1)%2) << ", " << ((share>>2)%2) << ((share>>3)%2) << std::endl;

      // Convert to format
      shares_x0[idx] = share%2;
      shares_x1[idx] = (share>>1)%2;
      shares_x0[idx + n] = (share>>2)%2;
      shares_x1[idx + n] = (share>>3)%2;
    }
  }

  // Execute conversion
  auto start = clock_start();
  store->heavy_ot(N, nbits, shares_x0, shares_x1, valid, bucket0, bucket1);
  std::cout << "heavy ot timing : " << sec_from(start) << std::endl;

  // recombine / test
  if (server_num == 0) {
    fmpz_t* bucket0_other; new_fmpz_array(&bucket0_other, nbits);
    fmpz_t* bucket1_other; new_fmpz_array(&bucket1_other, nbits);
    recv_fmpz_batch(serverfd, bucket0_other, nbits);
    recv_fmpz_batch(serverfd, bucket1_other, nbits);
    fmpz_t tmp; fmpz_init(tmp);
    for (unsigned int j = 0; j < nbits; j++) {
      fmpz_add(tmp, bucket0[j], bucket0_other[j]);
      fmpz_mod(tmp, tmp, Int_Modulus);
      std::cout << "bucket0[" << j << "] total = " << fmpz_get_ui(tmp) << std::endl;
      fmpz_add(tmp, bucket1[j], bucket1_other[j]);
      fmpz_mod(tmp, tmp, Int_Modulus);
      std::cout << "bucket1[" << j << "] total = " << fmpz_get_ui(tmp) << std::endl;
    }
    fmpz_clear(tmp);
    clear_fmpz_array(bucket0_other, nbits);
    clear_fmpz_array(bucket1_other, nbits);
  } else {
    send_fmpz_batch(serverfd, bucket0, nbits);
    send_fmpz_batch(serverfd, bucket1, nbits);
  }

  clear_fmpz_array(bucket0, nbits);
  clear_fmpz_array(bucket1, nbits);
  delete[] shares_x0;
  delete[] shares_x1;
  delete[] valid;
}


void test_cmp(
    const int server_num, const int serverfd, CorrelatedStore* store) {
  // setup
  const size_t N = 10;
  fmpz_t* val0; new_fmpz_array(&val0, N);
  fmpz_t* val1; new_fmpz_array(&val1, N);

  // "Very small" compared to modulus
  const size_t range = 1000;
  fmpz_t range_f; fmpz_init_set_si(range_f, range);
  for (unsigned int i = 0; i < N; i++) {
    fmpz_randm(val0[i], seed, range_f);
    fmpz_sub_si(val0[i], val0[i], range / 2);
    fmpz_mod(val0[i], val0[i], Int_Modulus);

    fmpz_randm(val1[i], seed, range_f);
    fmpz_sub_si(val1[i], val1[i], range / 2);
    fmpz_mod(val1[i], val1[i], Int_Modulus);
  }
  // fmpz_set_si(val0[0], -(100 + 10*server_num));
  // fmpz_mod(val0[0], val0[0], Int_Modulus);
  // fmpz_set_si(val1[0], (10 + server_num));
  // fmpz_mod(val1[0], val1[0], Int_Modulus);

  // run
  const bool* const larger = store->abs_cmp(N, val0, val1);

  // test
  if (server_num == 0) {
    fmpz_t* val0_other; new_fmpz_array(&val0_other, N);
    fmpz_t* val1_other; new_fmpz_array(&val1_other, N);
    recv_fmpz_batch(serverfd, val0_other, N);
    recv_fmpz_batch(serverfd, val1_other, N);

    fmpz_t v0; fmpz_init(v0);
    fmpz_t v1; fmpz_init(v1);
    fmpz_t half; fmpz_init(half);
    fmpz_cdiv_q_ui(half, Int_Modulus, 2);

    fmpz_t v0_tmp; fmpz_init(v0_tmp);
    fmpz_t v1_tmp; fmpz_init(v1_tmp);

    for (unsigned int i = 0; i < N; i++) {
      fmpz_add(v0, val0[i], val0_other[i]);
      fmpz_mod(v0, v0, Int_Modulus);
      fmpz_set(v0_tmp, v0);

      fmpz_add(v1, val1[i], val1_other[i]);
      fmpz_mod(v1, v1, Int_Modulus);
      fmpz_set(v1_tmp, v1);

      if (fmpz_cmp(v0_tmp, half) > 0) {  // > N/2, negate (abs)
        fmpz_sub(v0_tmp, Int_Modulus, v0);
      }
      if (fmpz_cmp(v1_tmp, half) > 0) {  // > N/2, negate (abs)
        fmpz_sub(v1_tmp, Int_Modulus, v1_tmp);
      }
      bool actual = fmpz_cmp(v0_tmp, v1_tmp) < 0;

      std::cout << i << ": |" << get_fsigned(v0, Int_Modulus) << (larger[i] ? "| < |" : "| > |") << get_fsigned(v1, Int_Modulus) << "|, \tactual: " << (actual ? "<" : ">") << std::endl;
      // std::cout << "(" << fmpz_get_ui(val0[i]) << " + " << fmpz_get_ui(val0_other[i]) << ") = " << get_fsigned(v0, Int_Modulus) << " vs " << get_fsigned(v1, Int_Modulus) << " = (" << fmpz_get_ui(val1[i]) << " + " << fmpz_get_ui(val1_other[i]) << ")\n";
    }
    fmpz_clear(v0);
    fmpz_clear(v1);
    fmpz_clear(half);
    clear_fmpz_array(val0_other, N);
    clear_fmpz_array(val1_other, N);
  } else {
    send_fmpz_batch(serverfd, val0, N);
    send_fmpz_batch(serverfd, val1, N);
  }

  clear_fmpz_array(val0, N);
  clear_fmpz_array(val1, N);
}

void test_cmp_bit(
    const int server_num, const int serverfd, CorrelatedStore* store) {
  // Setup values
  const size_t bits = 2;
  const size_t max = 1 << bits;
  const size_t N = max * max;
  fmpz_t* x_bits; new_fmpz_array(&x_bits, N * bits);
  fmpz_t* y_bits; new_fmpz_array(&y_bits, N * bits);
  if (server_num == 0) {
    fmpz_t* x_bits_other; new_fmpz_array(&x_bits_other, N * bits);
    fmpz_t* y_bits_other; new_fmpz_array(&y_bits_other, N * bits);
    for (unsigned int x = 0; x < max; x++) {
      for (unsigned int y = 0; y < max; y++) {
        // std::cout << "case x = " << x << ", y = " << y << std::endl;
        // if (x == y) continue;  // should be 0, skip for smaller logging for now
        for (unsigned int b = 0; b < bits; b++) {
          size_t idx = bits * (x * max + y) + b;
          fmpz_randm(x_bits[idx], seed, Int_Modulus);
          fmpz_sub(x_bits_other[idx], Int_Modulus, x_bits[idx]);
          fmpz_add_ui(x_bits_other[idx], x_bits_other[idx], (x >> b) % 2);

          fmpz_randm(y_bits[idx], seed, Int_Modulus);
          fmpz_sub(y_bits_other[idx], Int_Modulus, y_bits[idx]);
          fmpz_add_ui(y_bits_other[idx], y_bits_other[idx], (y >> b) % 2);

          // std::cout << " bit " << b << ", idx = " << idx << std::endl;
          // std::cout << "  x bit: " << ((x >> b) % 2) << " = " << fmpz_get_ui(x_bits[idx]) << " + " << fmpz_get_ui(x_bits_other[idx]) << std::endl;
          // std::cout << "  y bit: " << ((y >> b) % 2) << " = " << fmpz_get_ui(y_bits[idx]) << " + " << fmpz_get_ui(y_bits_other[idx]) << std::endl;
        }
      }
    }
    send_fmpz_batch(serverfd, x_bits_other, N * bits);
    send_fmpz_batch(serverfd, y_bits_other, N * bits);
    clear_fmpz_array(x_bits_other, N * bits);
    clear_fmpz_array(y_bits_other, N * bits);
  } else {
    recv_fmpz_batch(serverfd, x_bits, N * bits);
    recv_fmpz_batch(serverfd, y_bits, N * bits);
  }

  // Eval cmp
  fmpz_t* ans = store->cmp_bit(N, bits, x_bits, y_bits);

  // Check values
  if (server_num == 0) {
    fmpz_t* ans_other; new_fmpz_array(&ans_other, N);
    recv_fmpz_batch(serverfd, ans_other, N);

    fmpz_t tmp; fmpz_init(tmp);
    for (unsigned int x = 0; x < max; x++) {
      for (unsigned int y = 0; y < max; y++) {
        size_t idx = x * max + y;
        fmpz_add(tmp, ans[idx], ans_other[idx]);
        fmpz_mod(tmp, tmp, Int_Modulus);
        std::cout << "is " << x << " < " << y << "? " << (fmpz_is_one(tmp) ? "yes" : "no" ) << " (" << fmpz_get_ui(ans[idx]) << " + " << fmpz_get_ui(ans_other[idx]) << " = " << fmpz_get_ui(tmp) << ") " << std::endl;
      }
    }
    clear_fmpz_array(ans_other, N);
    fmpz_clear(tmp);
  } else {
    send_fmpz_batch(serverfd, ans, N);
  }


  clear_fmpz_array(ans, N);
  clear_fmpz_array(x_bits, N * bits);
  clear_fmpz_array(y_bits, N * bits);
}

void runServerTest(const int server_num, const int serverfd) {

  // init
  OT_Wrapper* ot0 = new OT_Wrapper(server_num == 0 ? nullptr : "127.0.0.1", 60051);
  OT_Wrapper* ot1 = new OT_Wrapper(server_num == 1 ? nullptr : "127.0.0.1", 60052);
  CorrelatedStore* store = new CorrelatedStore(serverfd, server_num, ot0, ot1, batch_size, lazy, do_fork);

  store->maybeUpdate();
  std::cout << std::endl;

  // random adjusting. different numbers adjust seed.
  // Can also do different per server
  fmpz_t tmp;
  fmpz_init(tmp);
  size_t rand_adjust;
  if (server_num == 0) {
    rand_adjust = 1;
  } else {
    rand_adjust = 5;
  }
  for (unsigned int i = 0; i < rand_adjust; i++)
    fmpz_randm(tmp, seed, Int_Modulus);
  fmpz_clear(tmp);

  // testHeavyOT(server_num, serverfd, store);

  // test_cmp(server_num, serverfd, store);

  test_cmp_bit(server_num, serverfd, store);
  
  delete ot0;
  delete ot1;
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

// 00, 01 = 0, 10 = 1, 11 = 0 (-1)
int main(int argc, char* argv[]) {
  init_constants();

  int server_num = -1;
  if(argc >= 2){
    server_num = atoi(argv[1]);
  }
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

  clear_constants();
}
