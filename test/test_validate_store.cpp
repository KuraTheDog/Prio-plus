#undef NDEBUG
#include <assert.h>
#include <iostream>

#include "utils_test_connect.h"
#include "../constants.h"
#include "../correlated_validated.h"
#include "../fmpz_utils.h"
#include "../net_share.h"
#include "../share.h"

#include "../interp.h"

const size_t rand_adjust = 5;
const bool do_fork = false;  // False for leaks debug
const bool lazy = true;

void multiServerRun(const int server_num, const int serverfd, const size_t N) {
  ValidateCorrelatedStore store(serverfd, server_num, N, lazy, do_fork);

  std::cout << "Setup " << server_num << std::endl;
  // Lazy pre-generate inputs
  if (server_num == 0) {
    // Gen
    DaBit* bits0[N];
    DaBit* bits1[N];
    AltTriple* trips0[N];
    AltTriple* trips1[N];
    for (unsigned int i = 0; i < N; i++) {
      bits0[i] = new DaBit();
      bits1[i] = new DaBit();
      makeLocalDaBit(bits0[i], bits1[i]);

      trips0[i] = new AltTriple();
      trips1[i] = new AltTriple();
      NewAltTriples(trips0[i], trips1[i]);
    }
    // Send
    send_DaBit_batch(serverfd, bits1, N);
    for (unsigned int i = 0; i < N; i++) {
      send_AltTriple(serverfd, trips1[i]);

      store.addUnverified(bits0[i], trips0[i]);

      delete bits1[i];
      delete trips1[i];
    }
  } else {
    DaBit* bits1[N];
    AltTriple* trips1[N];
    for (unsigned int i = 0; i < N; i++) {
      bits1[i] = new DaBit();
      trips1[i] = new AltTriple();
    }
    recv_DaBit_batch(serverfd, bits1, N);
    for (unsigned int i = 0; i < N; i++) {
       recv_AltTriple(serverfd, trips1[i]);
       store.addUnverified(bits1[i], trips1[i]);
    }
  }
  std::cout << "Run " << server_num << std::endl;

  store.batchValidate();

  std::cout << "Assert " << server_num << ", NumVerified = " << store.numVerified() << std::endl;

  assert(store.numVerified() == N);

}


void serverTest(const size_t N) {
  std::cout << "Running server test" << std::endl;
  int sockfd = init_receiver();

  pid_t pid = fork();
  if (pid == 0) {
    int cli_sockfd = init_sender();

    multiServerRun(0, cli_sockfd, N);

    close(cli_sockfd);
  } else {
    int newsockfd = accept_receiver(sockfd);

    fmpz_t tmp;
    fmpz_init(tmp);
    for (unsigned int i = 0; i < rand_adjust; i++)
      fmpz_randm(tmp, seed, Int_Modulus);
    fmpz_clear(tmp);

    multiServerRun(1, newsockfd, N);

    close(newsockfd);
  }

  close(sockfd);
}

int main(int argc, char* argv[]) {
  init_constants();


  serverTest(4);


  RootManager(1).clearCache();
  clear_constants();
}
