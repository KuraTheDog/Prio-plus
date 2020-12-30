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


/* Shares one server has for checking multiplcation gates
   See share_polynomial.
   Should satisfy f * g = h
   N = NextPowerOf2(M + 1), where M is the number of mul circuits.

   h_points is evaluated on 2N-th roots of unity, excluding points
   that are also N-th roots of unity. So all odd points.
*/
struct client_packet {
    size_t N;            // N = length of h_points = # mult gates.
    size_t NWires;       // Num wire shares = # mult gates + # input gates
    fmpz_t* WireShares;  // share of input/mulgate (output) wires
    fmpz_t f0_s;         // share of f(0) = pointsF[0] = u_0
    fmpz_t g0_s;         // share of g(0) = pointsG[0] = v_0
    fmpz_t h0_s;         // share of h(0)
    fmpz_t* h_points;    // h evaluated on odd 2N-th roots of unity

    ~client_packet() {
        fmpz_clear(f0_s);
        fmpz_clear(g0_s);
        fmpz_clear(h0_s);
        clear_fmpz_array(WireShares, NWires);
        clear_fmpz_array(h_points, N);
    }

    void print() {
        std::cout << " N = " << N << std::endl;
        std::cout << " NWires = " << NWires << std::endl;
        std::cout << " WireShares = {";
        for (int i = 0; i < NWires; i++) {
            if (i > 0)
                std::cout << ", ";
            fmpz_print(WireShares[i]);
        }
        std::cout << "}" << std::endl;
        std::cout << " f0_s = "; fmpz_print(f0_s); std::cout << std::endl;
        std::cout << " g0_s = "; fmpz_print(g0_s); std::cout << std::endl;
        std::cout << " h0_s = "; fmpz_print(h0_s); std::cout << std::endl;
        std::cout << " h_points = {";
        for (int i = 0; i < N; i++) {
            if (i > 0)
                std::cout << ", ";
            fmpz_print(h_points[i]);
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

void init_client_packet(ClientPacket &p, int N, int NumMulInpGates);

void SplitShare(fmpz_t val, fmpz_t A, fmpz_t B);

BeaverTriple* NewBeaverTriple();

BeaverTripleShare* BeaverTripleShares(BeaverTriple* inp);

#endif
