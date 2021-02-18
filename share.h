#ifndef SHARE_H
#define SHARE_H

#include <iostream>

// For fmpz types
#include "fmpz_utils.h"

struct Cor {
    fmpz_t D;
    fmpz_t E;

    Cor(){
        fmpz_init(D);
        fmpz_init(E);
    }

    ~Cor(){
        fmpz_clear(D);
        fmpz_clear(E);
    }
};

struct CorShare {
    fmpz_t shareD;
    fmpz_t shareE;

    CorShare(){
        fmpz_init(shareD);
        fmpz_init(shareE);
    }

    ~CorShare(){
        fmpz_clear(shareD);
        fmpz_clear(shareE);
    }
};

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

// I kind of want ot just combine this with BeaverTriple, since it's the same
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

struct BooleanBeaverTriple {
    bool a;
    bool b;
    bool c;

    BooleanBeaverTriple(){}

    BooleanBeaverTriple(const bool a, const bool b, const bool c)
    : a(a), b(b), c(c) {}
};

/* Shares one server has for checking multiplcation gates
   See share_polynomial.
   Should satisfy f * g = h
   N = NextPowerOf2(M + 1), where M is the number of mul circuits.

   h_points is evaluated on 2N-th roots of unity, excluding points
   that are also N-th roots of unity. So all odd points.
*/
struct ClientPacket {
    const size_t N;            // N = length of h_points = # mult gates.
    const size_t NWires;       // Num wire shares = # mult gates + # input gates
    fmpz_t* WireShares;  // share of input/mulgate (output) wires
    fmpz_t f0_s;         // share of f(0) = pointsF[0] = u_0
    fmpz_t g0_s;         // share of g(0) = pointsG[0] = v_0
    fmpz_t h0_s;         // share of h(0)
    fmpz_t* h_points;    // h evaluated on odd 2N-th roots of unity
    BeaverTripleShare* triple_share;

    ClientPacket(const size_t N, const size_t NumMulInpGates)
    : N(N), NWires(NumMulInpGates) {
        new_fmpz_array(&WireShares, NWires);
        fmpz_init(f0_s);
        fmpz_init(g0_s);
        fmpz_init(h0_s);
        new_fmpz_array(&h_points, N);
        triple_share = new BeaverTripleShare();
    }

    ~ClientPacket() {
        clear_fmpz_array(WireShares, NWires);
        fmpz_clear(f0_s);
        fmpz_clear(g0_s);
        fmpz_clear(h0_s);
        clear_fmpz_array(h_points, N);
        delete triple_share;
    }

    void print() const {
        std::cout << " N = " << N << std::endl;
        std::cout << " NWires = " << NWires << std::endl;
        std::cout << " WireShares = {";
        for (unsigned int i = 0; i < NWires; i++) {
            if (i > 0)
                std::cout << ", ";
            fmpz_print(WireShares[i]);
        }
        std::cout << "}" << std::endl;
        std::cout << " f0_s = "; fmpz_print(f0_s); std::cout << std::endl;
        std::cout << " g0_s = "; fmpz_print(g0_s); std::cout << std::endl;
        std::cout << " h0_s = "; fmpz_print(h0_s); std::cout << std::endl;
        std::cout << " h_points = {";
        for (unsigned int i = 0; i < N; i++) {
            if (i > 0)
                std::cout << ", ";
            fmpz_print(h_points[i]);
        }
        std::cout << "}" << std::endl;
    }
};

// Unused?
/*
struct ClientSubmission {
    fmpz_t* vals;
    BeaverTripleShare triple;
};
*/

struct DaBit {
    fmpz_t bp;   // [b]_p, mod p
    bool b2;     // [b]_2, bit

    DaBit() {
        fmpz_init(bp);
    }

    ~DaBit() {
        fmpz_clear(bp);
    }

    void print() const {
        std::cout << " [b]_p = "; fmpz_print(bp); std::cout << std::endl;
        std::cout << " [b]_2 = " << b2 << std::endl;
    }
};

struct EdaBit {
    fmpz_t r;  // [r]_p 
    const size_t n;
    bool* b;   // {[b]_2}

    EdaBit(const size_t n) : n(n) {
        fmpz_init(r);
        b = new bool[n];
    }

    ~EdaBit() {
        fmpz_clear(r);
        delete[] b;
    }

    // Get b as an int. For easier math / xor comparison.
    int get_int_b() const {
        int pow = 1, ans = 0;
        for (unsigned int i = 0; i < n; i++) {
            if (b[i])
                ans += pow;
            pow *= 2;
        }
        return ans;
    }

    void print() const {
        std::cout << " r = "; fmpz_print(r); std::cout << std::endl;
        fmpz_t x; fmpz_init(x);
        fmpz_from_bool_array(x, b, n);
        std::cout << " b = "; fmpz_print(x); std::cout << std::endl;
        fmpz_clear(x);
    }
};

void SplitShare(const fmpz_t val, fmpz_t A, fmpz_t B);

BeaverTriple* NewBeaverTriple();

void BeaverTripleShares(const BeaverTriple* const inp,
                        BeaverTripleShare* const out0,
                        BeaverTripleShare* const out1);

void makeLocalDaBit(DaBit* const bit0, DaBit* const bit1);
void makeLocalEdaBit(EdaBit* const ebit0, EdaBit* const ebit1, const size_t n);

#endif
