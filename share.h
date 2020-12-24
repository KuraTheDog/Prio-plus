#ifndef SHARE_H
#define SHARE_H

// For fmpz types
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

void init_client_packet(ClientPacket &p, int N);

void SplitShare(fmpz_t val, fmpz_t A, fmpz_t B);

BeaverTriple* NewBeaverTriple();

BeaverTripleShare* BeaverTripleShares(BeaverTriple* inp);

#endif