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

  // Build out of shares
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
  const size_t NMul;  // # mul gates = num wire shares
  const size_t N;     // N = length of h_points = # mult gates.
  fmpz_t* MulShares;   // share of mulgate (output) wires, aka h(even)
  fmpz_t f0_s;      // share of f(0) = pointsF[0] = u_0
  fmpz_t g0_s;      // share of g(0) = pointsG[0] = v_0
  fmpz_t h0_s;      // share of h(0)
  fmpz_t* h_points;   // h evaluated on odd 2N-th roots of unity
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

  void print() const;
};

struct DaBit {
  fmpz_t bp;   // [b]_p, mod p
  bool b2;   // [b]_2, bit

  DaBit() {
    fmpz_init(bp);
  }

  ~DaBit() {
    fmpz_clear(bp);
  }

  void print() const {
    std::cout << " [b]_p = " << fmpz_get_ui(bp) << "\n";
    std::cout << " [b]_2 = " << b2 << std::endl;
  }
};

// Random [r] such that [r] < M, for some `b`-bit modulus.
// Has shares [rB] of the bits of r.
// This version has arithmetic shares.
struct BitShare {
  const size_t b;
  fmpz_t r;
  fmpz_t* rB;  // Shares of the bits of r

  BitShare(const size_t b) : b(b) {
    fmpz_init(r);
    new_fmpz_array(&rB, b);
  }

  ~BitShare() {
    fmpz_clear(r);
    clear_fmpz_array(rB, b);
  }
};

// Deprecated. Multiple DaBits use less rounds
struct [[deprecated("Multiple Dabits requires less rounds")]] EdaBit {
  fmpz_t r;    // [r]_p
  const size_t n;  // number of bits, length of b
  bool* b;     // {[b]_2}

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

/*
With [num_inputs] each providing [num_values] values, accumulate into ans if valid.
values: size [num_inputs * num_values], [all values for input0, all for input1, ...]
valid: size [num_inputs], if false ignore that input
ans: size [num_values], ans[i] is modular sum of all input's value[i]
returns number of valid.
Here for convenience (to access mod), and since this is usually used to sum shares
*/
size_t accumulate(const size_t num_inputs, const size_t num_values,
                  const fmpz_t* const values, const bool* const valid,
                  fmpz_t* const ans);

// Note: A can't be val, but B can.
void SplitShare(const fmpz_t val, fmpz_t A, fmpz_t B);

void makeLocalBoolTriple(BooleanBeaverTriple* const out0, BooleanBeaverTriple* const out1);
void makeLocalTriple(BeaverTriple* const out0, BeaverTriple* const out1);
void makeLocalAltTriple(AltTriple* const out0, AltTriple* const out1);

void makeLocalDaBit(DaBit* const bit0, DaBit* const bit1);
// void makeLocalEdaBit(EdaBit* const ebit0, EdaBit* const ebit1, const size_t n);

#endif
