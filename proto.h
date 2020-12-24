/*
Oblivious Transfer code.
For BITSUM and INTSUM.
*/

#ifndef PROTO_H
#define PROTO_H

#include <emp-ot/emp-ot.h>
#include <emp-tool/emp-tool.h>
#include <iostream>
using namespace emp;

uint64_t bitsum_ot_sender(NetIO *io,bool *shares, bool *valid, int n);

uint64_t intsum_ot_sender(NetIO *io,uint32_t *shares, bool *valid, int n, int num_bits);

uint64_t bitsum_ot_receiver(NetIO *io,bool *shares, int n);

uint64_t intsum_ot_receiver(NetIO *io, uint32_t *shares, int n, int num_bits);

#endif
