#ifndef SHARE_H
#define SHARE_H

#include <vector>
#include <string>
#include <iostream>
#include <gmpxx.h>

// #include "circuit.h"
#include "fmpz_utils.h"


struct Cor {
    fmpz_t D;
    fmpz_t E;
};

struct CorShare {
    fmpz_t shareD;
    fmpz_t shareE;
};

struct client_packet {
    size_t N;  // length of WireShares and h_points
    fmpz_t* WireShares;
    fmpz_t f0_s;
    fmpz_t g0_s;
    fmpz_t h0_s;
    fmpz_t* h_points;
};

typedef struct client_packet* ClientPacket;

struct BeaverTriple {
    fmpz_t A;
    fmpz_t B;
    fmpz_t C;

    BeaverTriple(){
        fmpz_init(A);
        fmpz_init(B);
        fmpz_init(C);
    }

    ~BeaverTriple(){
        fmpz_clear(A);
        fmpz_clear(B);
        fmpz_clear(C);
    }
};

struct BeaverTripleShare {
    fmpz_t shareA;
    fmpz_t shareB;
    fmpz_t shareC;

    BeaverTripleShare(){
        fmpz_init(shareA);
        fmpz_init(shareB);
        fmpz_init(shareC);
    }

    ~BeaverTripleShare(){
        fmpz_clear(shareA);
        fmpz_clear(shareB);
        fmpz_clear(shareC);
    }
};

struct ClientSubmission {
    fmpz_t* vals;
    BeaverTripleShare triple;
};

fmpz_t Int_Modulus;
flint_rand_t seed;

fmpz_t* SplitShare(fmpz_t val) {
    fmpz_t* ans = (fmpz_t*) malloc(2*sizeof(fmpz_t));

    fmpz_init(ans[0]);
    fmpz_init(ans[1]);

    fmpz_randm(ans[0],seed,Int_Modulus);
    fmpz_sub(ans[1],val,ans[0]);
    fmpz_mod(ans[1],ans[1],Int_Modulus);

    return ans;
}

void init_client_packet(ClientPacket &p, int N){
    p = new client_packet();

    p->N = N;

    new_fmpz_array(&p->WireShares, N);
    
    fmpz_init(p->f0_s);
    fmpz_init(p->g0_s);
    fmpz_init(p->h0_s);

    new_fmpz_array(&p->h_points, N);
    
}

void SplitShare(fmpz_t val, fmpz_t A, fmpz_t B){
    fmpz_randm(A,seed,Int_Modulus);
    fmpz_sub(B,val,A);
    fmpz_mod(B,B,Int_Modulus);
}

BeaverTriple* NewBeaverTriple() {
    BeaverTriple* out = new BeaverTriple();
    fmpz_randm(out->A,seed,Int_Modulus);
    fmpz_randm(out->B,seed,Int_Modulus);
    fmpz_mul(out->C,out->A,out->B);
    fmpz_mod(out->C,out->C,Int_Modulus);
    return out;
}

BeaverTripleShare* BeaverTripleShares(BeaverTriple* inp) {
    BeaverTripleShare* out = new BeaverTripleShare[2];

    SplitShare(inp->A,out[0].shareA,out[1].shareA);
    SplitShare(inp->B,out[0].shareB,out[1].shareB);
    SplitShare(inp->C,out[0].shareC,out[1].shareC);

    return out;
}

#endif