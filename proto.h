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

BooleanBeaverTriple* gen_boolean_beaver_triples(const int server_num, const unsigned int m, NetIO* const io0, NetIO* const io1);

BeaverTriple* generate_beaver_triple(const int serverfd, const int server_num, NetIO* const io0, NetIO* const io1);

BeaverTriple* generate_beaver_triple_lazy(const int serverfd, const int server_num);

#endif
