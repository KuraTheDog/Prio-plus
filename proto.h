/**
Oblivious Transfer code.
For BITSUM and INTSUM.
*/

#ifndef PROTO_H
#define PROTO_H

#include <emp-ot/emp-ot.h>
#include <emp-tool/emp-tool.h>
#include <iostream>

#include "share.h"

using namespace emp;

uint64_t bitsum_ot_sender(NetIO* const io, const bool* const shares, const bool* const valid, const size_t n);

uint64_t intsum_ot_sender(NetIO* const io, const uint64_t* const shares, const bool* const valid, const size_t n, const size_t num_bits);

uint64_t bitsum_ot_receiver(NetIO* const io, const bool* const shares, const size_t n);

uint64_t intsum_ot_receiver(NetIO* const io, const uint64_t* const shares, const size_t n, const size_t num_bits);

BooleanBeaverTriple* gen_boolean_beaver_triples(const int server_num, const int m, NetIO* const io0, NetIO* const io1);

BeaverTriple* generate_beaver_triple(const int serverfd, const int server_num, NetIO* const io0, NetIO* const io1);

struct SHEkeys {
  // const size_t pk_length;
  const int server_num;

  const int my_pk;
  const int my_sk;
  int other_pk;

  SHEkeys(const int server_num) 
  : server_num(server_num)
  , my_pk(server_num == 0 ? 888800 : 888811)
  , my_sk(server_num == 0 ? 777700 : 777711)
  {
    // std::cout << server_num << " made shekeys" << std::endl;
    // std::cout << server_num << " my pk: " << my_pk << std::endl;
    // std::cout << server_num << " my sk: " << my_sk << std::endl;
  }

  ~SHEkeys() {}

  void exchange_pk(const int serverfd);
};

BeaverTriple* generate_beaver_triple_she(const int serverfd, const int server_num, const SHEkeys* keys);

#endif
