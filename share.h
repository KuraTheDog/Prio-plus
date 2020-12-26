#ifndef SHARE_H
#define SHARE_H

#include <iostream>

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

    void print() {
        std::cout << " N = " << N << std::endl;
        std::cout << " WireShares = {";
        for (int i = 0; i < N; i++) {
            if (i > 0)
                std::cout << ", ";
            std::cout << fmpz_print(WireShares[i]);
        }
        std::cout << "}" << std::endl;
        std::cout << " f0_s = " << f0_s << std::endl;
        std::cout << " g0_s = " << g0_s << std::endl;
        std::cout << " h0_s = " << h0_s << std::endl;
        std::cout << " h_points = {";
        for (int i = 0; i < N; i++) {
            if (i > 0)
                std::cout << ", ";
            std::cout << fmpz_print(h_points[i]);
        }
        std::cout << "}" << std::endl;
    }
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
