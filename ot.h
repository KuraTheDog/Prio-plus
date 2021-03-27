/*
Oblivious Transfer code
For share conversion
*/

#ifndef PROTO_H
#define PROTO_H

#include <emp-ot/emp-ot.h>
#include <emp-tool/emp-tool.h>
#include <iostream>
#include <queue>

#include "share.h"

#define EMP_IKNP 1

#define OT_TYPE EMP_IKNP

struct OT_Wrapper {
  emp::NetIO* const io;
  emp::IKNP<emp::NetIO>* const ot;

  OT_Wrapper(const char* address, const int port);
  ~OT_Wrapper();

  void send(const uint64_t* const data0, const uint64_t* const data1,
            const size_t length);
  void recv(uint64_t* const data, const bool* b, const size_t length);
};

uint64_t bitsum_ot_sender(OT_Wrapper* const ot, const bool* const shares, const bool* const valid, const size_t n);
uint64_t bitsum_ot_receiver(OT_Wrapper* const ot, const bool* const shares, const size_t n);

// Batched version, for multiple values
// shares: #shares * #values, as (s0v0, s0v1, s0v2, s1v0, ...)
// valid: validity of share i
// bits: length of value j
uint64_t* intsum_ot_sender(OT_Wrapper* const ot, const uint64_t* const shares,
                           const bool* const valid, const size_t* const num_bits,
                           const size_t num_shares, const size_t num_values);
uint64_t* intsum_ot_receiver(OT_Wrapper* const ot, const uint64_t* const shares,
                             const size_t* const num_bits,
                             const size_t num_shares, const size_t num_values);

std::queue<BooleanBeaverTriple*> gen_boolean_beaver_triples(const int server_num, const unsigned int m, OT_Wrapper* const ot0, OT_Wrapper* const ot1);

BeaverTriple* generate_beaver_triple(const int serverfd, const int server_num, OT_Wrapper* const ot0, OT_Wrapper* const ot1);

BeaverTriple* generate_beaver_triple_lazy(const int serverfd, const int server_num);

#endif
