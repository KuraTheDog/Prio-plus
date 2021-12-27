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
void runServerTest(const int server_num, const int serverfd) {

  // set up
  OT_Wrapper* ot0 = new OT_Wrapper(server_num == 0 ? nullptr : "127.0.0.1", 60051);
  OT_Wrapper* ot1 = new OT_Wrapper(server_num == 1 ? nullptr : "127.0.0.1", 60052);
  CorrelatedStore* store = new CorrelatedStore(serverfd, server_num, ot0, ot1, batch_size, lazy, do_fork);

  store->maybeUpdate();
  std::cout << std::endl;

  const size_t nbits = 4;
  const size_t N = 8;
  const size_t n = nbits * N;

  bool* const shares_x0 = new bool[2 * n];
  bool* const shares_x1 = new bool[2 * n];
  fmpz_t* bucket0; new_fmpz_array(&bucket0, nbits);
  fmpz_t* bucket1; new_fmpz_array(&bucket1, nbits);
  bool* const valid = new bool[N]; memset(valid, 1, N);

  // Set up values
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

      // Convert to format. Is this right?
      shares_x0[idx] = share%2;
      shares_x1[idx] = (share>>1)%2;
      shares_x0[idx + n] = (share>>2)%2;
      shares_x1[idx + n] = (share>>3)%2;
      // std::cout << " set shares_x0[" << idx << "] = " << shares_x0[idx] << std::endl;
      // std::cout << " set shares_x0[" << idx+n << "] = " << shares_x0[idx+n] << std::endl;
    }
  }

  // std::cout << "outside shares_x0: ";
  // for (unsigned int i = 0; i < 2 * n; i++) {
  //   std::cout << shares_x0[i];
  //   if (i == n - 1) std::cout << " ";
  // }
  // std::cout << std::endl;

  auto start = clock_start();
  store->heavy_ot(N, nbits, shares_x0, shares_x1, valid, bucket0, bucket1);
  std::cout << "heavy ot timing : " << sec_from(start) << std::endl;

  // recombine
  if (server_num == 0) {
    fmpz_t* bucket0_other; new_fmpz_array(&bucket0_other, nbits);
    fmpz_t* bucket1_other; new_fmpz_array(&bucket1_other, nbits);
    recv_fmpz_batch(serverfd, bucket0_other, nbits);
    recv_fmpz_batch(serverfd, bucket1_other, nbits);
    fmpz_t tmp; fmpz_init(tmp);
    for (unsigned int i = 0; i < nbits; i++) {
      fmpz_add(tmp, bucket0[i], bucket0_other[i]);
      fmpz_mod(tmp, tmp, Int_Modulus);
      std::cout << "bucket0[" << i << "] total = " << fmpz_get_ui(tmp) << std::endl;
      fmpz_add(tmp, bucket1[i], bucket1_other[i]);
      fmpz_mod(tmp, tmp, Int_Modulus);
      std::cout << "bucket1[" << i << "] total = " << fmpz_get_ui(tmp) << std::endl;
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
  clear_fmpz_array(bucket0, nbits);
  clear_fmpz_array(bucket1, nbits);
  delete[] valid;
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
