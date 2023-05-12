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
const bool lazy = true;

void setup(const int server_num, const int serverfd, const size_t N,
           ValidateCorrelatedStore& store) {
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

      store.addUnvalidated(bits0[i], trips0[i]);

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
       store.addUnvalidated(bits1[i], trips1[i]);
    }
  }
}

void test_altMult(const int server_num, const int serverfd,
                  const size_t N_make, const size_t N) {
  OT_Wrapper* ot0 = new OT_Wrapper(server_num == 0 ? nullptr : "127.0.0.1", 60051);
  OT_Wrapper* ot1 = new OT_Wrapper(server_num == 1 ? nullptr : "127.0.0.1", 60052);

  ValidateCorrelatedStore store(serverfd, server_num, ot0, ot1, N_make, lazy);

  // Setup: Not necessary, use of check_AltTriples should auto-gen.
  // setup(server_num, serverfd, N, store);
  store.printSizes();

  fmpz_t* x; new_fmpz_array(&x, N);
  bool use_validated[N];
  memset(use_validated, false, N);
  for (unsigned int i = 0; i < N; i++) {
    // 1 * 10, 2 * 20, ...
    fmpz_set_si(x[i], (i+1) * (1 + server_num * 9));
  }

  fmpz_t* z; new_fmpz_array(&z, N);
  store.multiplyAltShares(N, x, z, use_validated);
  swap_fmpz_batch(serverfd, z, N);
  for (unsigned int i = 0; i < N; i++) {
    assert(fmpz_equal_ui(z[i], (i+1)*(i+1) * 10));
  }

  clear_fmpz_array(z, N);
  delete ot0;
  delete ot1;
}

// TODO: also check checkDaBits
void test_batchValidate(const int server_num, const int serverfd,
                        const size_t N_make, const size_t N) {
  OT_Wrapper* ot0 = new OT_Wrapper(server_num == 0 ? nullptr : "127.0.0.1", 60051);
  OT_Wrapper* ot1 = new OT_Wrapper(server_num == 1 ? nullptr : "127.0.0.1", 60052);

  if (server_num == 0) std::cout << "Testing batch validate size: " << N << std::endl;
  ValidateCorrelatedStore store(serverfd, server_num, ot0, ot1, N_make, lazy);

  setup(server_num, serverfd, N, store);
  // store.printSizes();

  store.batchValidate();

  std::cout << "Assert " << server_num << ", Numvalidated = " << store.num_validated_dabits() << std::endl;

  assert(store.num_validated_dabits() == N);

  store.printSizes();

  delete ot0;
  delete ot1;
}


void serverTest(const size_t N) {
  std::cout << "Running server test" << std::endl;
  int sockfd = init_receiver();

  pid_t pid = fork();
  if (pid == 0) {
    int cli_sockfd = init_sender();

    test_altMult(0, cli_sockfd, N, N);

    test_batchValidate(0, cli_sockfd, N, N);

    close(cli_sockfd);
  } else {
    int newsockfd = accept_receiver(sockfd);

    fmpz_t tmp;
    fmpz_init(tmp);
    for (unsigned int i = 0; i < rand_adjust; i++)
      fmpz_randm(tmp, seed, Int_Modulus);
    fmpz_clear(tmp);

    test_altMult(1, newsockfd, N, N);

    test_batchValidate(1, newsockfd, N, N);

    close(newsockfd);
  }

  close(sockfd);
}

int main(int argc, char* argv[]) {
  init_constants();

  int N = 4;
  if (argc >= 2) {
    N = atoi(argv[1]);
  }

  serverTest(N);


  RootManager(1).clearCache();
  clear_constants();
}
