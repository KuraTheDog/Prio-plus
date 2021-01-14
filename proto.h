/**
Oblivious Transfer code.
For BITSUM and INTSUM.
*/

#ifndef PROTO_H
#define PROTO_H

#include <emp-ot/emp-ot.h>
#include <emp-tool/emp-tool.h>
#include <iostream>
using namespace emp;

// Boolean beaver triple
struct bbt {
    uint32_t x, y, z;

    bbt(bool xx, bool yy, bool zz){
        x = xx ? 1 : 0;
        y = yy ? 1 : 0;
        z = zz ? 1 : 0;
    }
};

uint64_t bitsum_ot_sender(NetIO *io,bool *shares, bool *valid, int n);

uint64_t intsum_ot_sender(NetIO *io,uint32_t *shares, bool *valid, int n, int num_bits);

uint64_t bitsum_ot_receiver(NetIO *io,bool *shares, int n);

uint64_t intsum_ot_receiver(NetIO *io, uint32_t *shares, int n, int num_bits);

vector<bbt> gen_boolean_beaver_triples(int server_num, int m);

#endif
