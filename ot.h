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

using namespace emp;

// Persistent OT with wrapper
struct OT_Wrapper {
  NetIO* const io;
  IKNP<NetIO>* const ot;

  OT_Wrapper(const char* address, const int port)
  : io(new NetIO(address, port, true))
  , ot(new IKNP<NetIO>(io))
  {}

  ~OT_Wrapper() {
    delete ot;
    delete io;
  }

  void send(const block* data0, const block* data1, const size_t length) {
    io->sync();
    ot->send(data0, data1, length);
    io->flush();
  }

  void recv(block* data, const bool* b, const size_t length) {
    io->sync();
    ot->recv(data, b, length);
    io->flush();
  }
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
