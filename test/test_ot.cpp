#undef NDEBUG
#include <assert.h>
#include <iostream>

#include "utils_test_connect.h"
#include "../constants.h"
#include "../net_share.h"
#include "../ot.h"

#define SERVER0_IP "127.0.0.1"
#define SERVER1_IP "127.0.0.1"

#define MOD 100
// #define MOD 0  // 2^64, default

// Config
const int N = 2;
const bool test_base = true;
const bool test_bit = true;
const bool test_int = true;

void run_server0(const size_t m) {
  const int server_num = 0;

  int cli_sockfd = init_sender();

  OT_Wrapper* ot0 = new OT_Wrapper(server_num == 0 ? nullptr : SERVER0_IP, 60051);
  OT_Wrapper* ot1 = new OT_Wrapper(server_num == 1 ? nullptr : SERVER1_IP, 60052);

  // General act like one input, all valid
  const bool valid[1] = {true};
  uint64_t b;

  if (test_base) {
    std::cout << "Testing basic send/recv" << std::endl;
    uint64_t* data0 = new uint64_t[N]; memset(data0, 0, N * sizeof(uint64_t));
    uint64_t* data0_1 = new uint64_t[N]; memset(data0_1, 0, N * sizeof(uint64_t));
    uint64_t* data1 = new uint64_t[N]; memset(data1, 0, N * sizeof(uint64_t));
    uint64_t* data1_1 = new uint64_t[N]; memset(data1_1, 0, N * sizeof(uint64_t));
    data0[0] = 10; data0[1] = 20;
    data0_1[0] = 11; data0_1[1] = 21;
    data1[0] = 1000; data1[1] = 2000;
    data1_1[0] = 1001; data1_1[1] = 2001;
    ot0->send(data0, data1, N, data0_1, data1_1);
  }

  if (test_bit) {
    std::cout << "Testing bit OT convert" << std::endl;
    const bool bitshares[1] = {(m / 2) % 2 ? true : false};
    // std::cout << "bit share0[0] = " << bitshares[0] << std::endl;
    const uint64_t a = bitsum_ot_sender(ot0, bitshares, valid, 1, MOD);
    recv_uint64(cli_sockfd, b);
    // std::cout << "bit have a = " << a] << std::endl;
    // std::cout << "bit got  b = " << b << std::endl;
    assert((a + b) % MOD == (m / 2) % 2);
  }

  if (test_int) {
    std::cout << "Testing int OT convert" << std::endl;
    // const uint64_t intshares[1][1] = {{m ^ 9}};
    uint64_t** const intshares = new uint64_t*[1];
    intshares[0] = new uint64_t[1];
    intshares[0][0] = m ^ 9;
    const size_t sizes[1] = {4};
    // std::cout << "int share0[0] = " << intshares[0][0] << std::endl;
    const uint64_t* const * const int_a = intsum_ot_sender(ot0,
        intshares, valid, sizes, 1, 1, MOD);
    delete[] intshares[0];
    delete[] intshares;
    recv_uint64(cli_sockfd, b);
    // std::cout << "int have a = " << int_a[0][0] << std::endl;
    // std::cout << "int got  b = " << b << std::endl;
    assert((int_a[0][0] + b) % MOD == m);
    delete[] int_a[0];
    delete[] int_a;
  }
  // cleanup

  delete ot0;
  delete ot1;
  close(cli_sockfd);
}

void run_server1(const size_t m) {
  const int server_num = 1;

  int sockfd = init_receiver();
  int newsockfd = accept_receiver(sockfd);

  OT_Wrapper* ot0 = new OT_Wrapper(server_num == 0 ? nullptr : SERVER0_IP, 60051);
  OT_Wrapper* ot1 = new OT_Wrapper(server_num == 1 ? nullptr : SERVER1_IP, 60052);

  if (test_base) {
    uint64_t* data = new uint64_t[N];
    uint64_t* data_1 = new uint64_t[N];
    bool* c = new bool[N]; memset(c, 0, N); c[0] = false; c[1] = true;
    ot0->recv(data, c, N, data_1);
    assert(data[0] == 10); assert(data_1[0] == 11);
    assert(data[1] == 2000); assert(data_1[1] == 2001);
    delete[] data;
    delete[] data_1;
  }

  if (test_bit) {
    const bool bitshares[1] = {m % 2 ? true : false};
    std::cout << "bit share1[0] = " << bitshares[0] << std::endl;
    const uint64_t b = bitsum_ot_receiver(ot0, bitshares, 1, MOD);
    send_uint64(newsockfd, b);
  }

  if (test_int) {
    uint64_t** const intshares = new uint64_t*[9];
    intshares[0] = new uint64_t[1];
    intshares[0][0] = 9;
    const size_t sizes[1] = {4};
    std::cout << "int share1[0] = " << intshares[0][0] << std::endl;
    const uint64_t* const * const int_b = intsum_ot_receiver(ot0,
        intshares, sizes, 1, 1, MOD);
    delete[] intshares[0];
    delete[] intshares;
    send_uint64(newsockfd, int_b[0][0]);
    delete[] int_b[0];
    delete[] int_b;
  }

  delete ot0;
  delete ot1;
  close(sockfd);
  close(newsockfd);
}

int main(int argc, char** argv){
  int m = 4;
  if (argc >= 2){
    m = atoi(argv[1]);
  }
  std::cout << "Input: " << m << "\n";

  init_constants();

  std::thread t0([m](){run_server0(m);});
  std::thread t1([m](){run_server1(m);});
  t0.join();
  t1.join();

  clear_constants();

  return 0;
}
