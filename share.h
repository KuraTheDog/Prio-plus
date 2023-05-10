#ifndef SHARE_H
#define SHARE_H

#include <iostream>

// For fmpz types
#include "fmpz_utils.h"

// For checker's triples
struct Cor {
    fmpz_t D;
    fmpz_t E;

    Cor(){
        fmpz_init(D);
        fmpz_init(E);
    }

    Cor(const Cor* const x, const Cor* const y);

    ~Cor(){
        fmpz_clear(D);
        fmpz_clear(E);
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

struct BooleanBeaverTriple {
    bool a;
    bool b;
    bool c;

    BooleanBeaverTriple(){}

    BooleanBeaverTriple(const bool a, const bool b, const bool c)
    : a(a), b(b), c(c) {}
};

// S0 has A, S1 has B, both have [C]
struct AltTriple {
    fmpz_t AB;
    fmpz_t C;
    AltTriple(){
        fmpz_init(AB);
        fmpz_init(C);
    }
    ~AltTriple(){
        fmpz_clear(AB);
        fmpz_clear(C);
    }
};

// Not inclusive, so next(2) = 4
unsigned int NextPowerOfTwo(const unsigned int n);

/* Shares one server has for checking multiplcation gates
   See share_polynomial.
   Should satisfy f * g = h
   N = NextPowerOf2(NMul), where M is the number of mul circuits.

   h_points is evaluated on 2N-th roots of unity, excluding points
   that are also N-th roots of unity. So all odd points.
*/
struct ClientPacket {
    const size_t NMul;    // # mul gates = num wire shares
    const size_t N;       // N = length of h_points = # mult gates.
    fmpz_t* MulShares;   // share of mulgate (output) wires, aka h(even)
    fmpz_t f0_s;          // share of f(0) = pointsF[0] = u_0
    fmpz_t g0_s;          // share of g(0) = pointsG[0] = v_0
    fmpz_t h0_s;          // share of h(0)
    fmpz_t* h_points;     // h evaluated on odd 2N-th roots of unity
    BeaverTriple* triple;

    ClientPacket(const size_t NMul)
    : NMul(NMul), N(NextPowerOfTwo(NMul)) {
        new_fmpz_array(&MulShares, NMul);
        fmpz_init(f0_s);
        fmpz_init(g0_s);
        fmpz_init(h0_s);
        new_fmpz_array(&h_points, N);
        triple = new BeaverTriple();
    }

    ~ClientPacket() {
        clear_fmpz_array(MulShares, NMul);
        fmpz_clear(f0_s);
        fmpz_clear(g0_s);
        fmpz_clear(h0_s);
        clear_fmpz_array(h_points, N);
        delete triple;
    }

    void print() const {
        std::cout << " NMul = " << NMul << "\n";
        std::cout << " N = " << N << "\n";
        std::cout << " MulShares = {";
        for (unsigned int i = 0; i < NMul; i++) {
            if (i > 0)
                std::cout << ", ";
            fmpz_print(MulShares[i]);
        }
        std::cout << "}\n";
        std::cout << " f0_s = "; fmpz_print(f0_s); std::cout << "\n";
        std::cout << " g0_s = "; fmpz_print(g0_s); std::cout << "\n";
        std::cout << " h0_s = "; fmpz_print(h0_s); std::cout << "\n";
        std::cout << " h_points = {";
        for (unsigned int i = 0; i < N; i++) {
            if (i > 0)
                std::cout << ", ";
            fmpz_print(h_points[i]);
        }
        std::cout << "}" << std::endl;
    }
};

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
        std::cout << " [b]_p = "; fmpz_print(bp); std::cout << "\n";
        std::cout << " [b]_2 = " << b2 << std::endl;
    }
};

// Deprecated. Multiple DaBits use less rounds
struct [[deprecated("Multiple Dabits requires less rounds")]] EdaBit {
    fmpz_t r;        // [r]_p
    const size_t n;  // number of bits, length of b
    bool* b;         // {[b]_2}

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
        std::cout << " r = " << fmpz_get_ui(r) << "\n";
        fmpz_t x; fmpz_init(x);
        fmpz_from_bool_array(x, b, n);
        std::cout << " b = "  << fmpz_get_ui(x) << std::endl;
        fmpz_clear(x);
    }
};

// Note: A can't be val, but B can.
void SplitShare(const fmpz_t val, fmpz_t A, fmpz_t B);

void NewBeaverTriples(BeaverTriple* const out0, BeaverTriple* const out1);
void NewAltTriples(AltTriple* const out0, AltTriple* const out1);

void makeLocalDaBit(DaBit* const bit0, DaBit* const bit1);
[[deprecated]] void makeLocalEdaBit(EdaBit* const ebit0, EdaBit* const ebit1, const size_t n);

#endif
