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

  // Communication costs (computed via monitoring network traffic)
  // Is bytes sent, since receiver also sends.
  // There are also one-time-setup costs, which are on server startup so not part of communication
  const size_t bytes_sender_start = 1;
  const size_t bytes_sender_per = 32;
  const size_t bytes_recver_per_block = 2048;
  const size_t bytes_recver_block_size = 128;

  OT_Wrapper(const char* const address, const int port);
  ~OT_Wrapper();

  // OT send/recv 64 bit data = datab selected from (data0, data) using b
  // If dataX_1 are set, is additional 64 bits added on, for 128 total.
  // I.e. (data, data_1) = (datab, datab_1) for choice b
  int send(const uint64_t* const data0, const uint64_t* const data1,
           const size_t length,
           const uint64_t* const data0_1 = nullptr, const uint64_t* const data1_1 = nullptr);
  int recv(uint64_t* const data, const bool* b, const size_t length,
           uint64_t* const data_1 = nullptr);
};

// mod 0 = default 2^64
uint64_t bitsum_ot_sender(
    OT_Wrapper* const ot, const bool* const shares, const bool* const valid,
    const size_t n, const size_t mod = 0);
uint64_t bitsum_ot_receiver(
    OT_Wrapper* const ot, const bool* const shares,
    const size_t n, const size_t mod = 0);

// Batched version, for multiple values per share
// shares and ret: matrix of shares x values, as [s0v0, s0v1, s0v2], [s1v0, ...]
// valid: validity of share i
// bits: length of value j
// Does not accumulate
const uint64_t* const * const intsum_ot_sender(
    OT_Wrapper* const ot, const uint64_t* const * const shares, const bool* const valid,
    const size_t* const num_bits, const size_t num_shares, const size_t num_values,
    const size_t mod = 0);
const uint64_t* const * const intsum_ot_receiver(
    OT_Wrapper* const ot, const uint64_t* const * const shares,
    const size_t* const num_bits, const size_t num_shares, const size_t num_values,
    const size_t mod = 0);

std::queue<const BooleanBeaverTriple* const> gen_boolean_beaver_triples(
    const int server_num, const unsigned int m,
    OT_Wrapper* const ot0, OT_Wrapper* const ot1);

const BeaverTriple* const generate_beaver_triple(
    const int serverfd, const int server_num,
    OT_Wrapper* const ot0, OT_Wrapper* const ot1);

const BeaverTriple* const generate_beaver_triple_lazy(
    const int serverfd, const int server_num);

#endif
