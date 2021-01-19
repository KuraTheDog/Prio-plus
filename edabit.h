#ifndef EDABIT_H
#define EDABIT_H

#include <iostream>

#include "constants.h"
#include "fmpz_utils.h"


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
    std::cout << " [b]_p = "; fmpz_print(bp); std::cout << std::endl;
    std::cout << " [b]_2 = " << b2 << std::endl;
  }
};

void makeLocalDaBit(DaBit* bit0, DaBit* bit1) {
  fmpz_t two;  // TODO: Is there a nicer way to do this?
  fmpz_init_set_si(two, 2);

  // random bit b
  fmpz_t bit;
  fmpz_init(bit);
  fmpz_randm(bit, seed, two);

  // random r
  fmpz_randm(bit0->bp, seed, Int_Modulus);

  // b - r
  fmpz_sub(bit1->bp, bit, bit0->bp);
  fmpz_mod(bit1->bp, bit1->bp, Int_Modulus);

  // r mod 2. true if odd.
  bit0->b2 = fmpz_is_odd(bit0->bp);

  // 1 xor (b - r) mod 2 = (1 + b - r) mod 2, true if (b-r) even
  bit1->b2 = fmpz_is_even(bit1->bp);

  fmpz_clear(two);
  fmpz_clear(bit);
}


/*
Single player's share
*/

struct EdaBit {
  fmpz_t r;  // [r]_p 
  const size_t n;
  bool* b;

  EdaBit(const size_t n) : n(n) {
    fmpz_init(r);
    b = new bool[n];
  }

  ~EdaBit() {
    fmpz_clear(r);
    delete[] b;
  }

  int get_int_b() const {
    int pow = 1, ans = 0;
    for (int i = 0; i < n; i++) {
      if (b[i])
        ans += pow;
      pow *= 2;
    }
    return ans;
  }

  void print() const {
    std::cout << " r = "; fmpz_print(r); std::cout << std::endl;
    std::cout << " b = " << get_int_b() << std::endl;
  }
};

void makeLocalEdaBit(EdaBit* ebit0, EdaBit* ebit1, const size_t n) {
  DaBit* bit0 = new DaBit();
  DaBit* bit1 = new DaBit();

  fmpz_zero(ebit0->r);
  fmpz_zero(ebit1->r);

  fmpz_t pow;
  fmpz_init_set_si(pow, 1);  // 2^i

  for (int i = 0; i < n; i++) {
    makeLocalDaBit(bit0, bit1);

    // std::cout << "bit0[" << i << "]\n";
    // bit0->print();
    // std::cout << "bit1[" << i << "]\n";
    // bit1->print();

    // if (bit0->b2)
    //   fmpz_add(ebit0->b, ebit0->b, pow);
    // if (bit1->b2)
    //   fmpz_add(ebit1->b, ebit1->b, pow);
    ebit0->b[i] = bit0->b2;
    ebit1->b[i] = bit1->b2;

    fmpz_addmul(ebit0->r, pow, bit0->bp);
    fmpz_mod(ebit0->r, ebit0->r, Int_Modulus);
    fmpz_addmul(ebit1->r, pow, bit1->bp);
    fmpz_mod(ebit1->r, ebit1->r, Int_Modulus);

    fmpz_mul_si(pow, pow, 2);
    fmpz_mod(pow, pow, Int_Modulus);
  }
  fmpz_clear(pow);
  delete bit0;
  delete bit1;
}

void makeEdaBit(EdaBit* ebit0, EdaBit* ebit1, const size_t n) {

}

// Todo: deal with beaver triples
// Get [x] + [y] = [z] + [carry], arrays of size n.
bool addBitShares(const bool* x, const bool* y, bool* z, const size_t n, const int server_num) {
  bool carry = false;
  bool x2, y2, prod;

  bool d, e, d_all, e_all;

  bool a = 0, b = 0, c = 0;  // Beaver triples. Todo: get from params.

  for (int i = 0; i < n; i++) {
    z[i] = carry ^ x[i] ^ y[i];
    x2 = carry ^ x[i];
    y2 = carry ^ y[i];

    d = x2 ^ a;
    e = y2 ^ b;
    // broadcast d and e. Get other d and e. Depends on otehr x, y
    d_all = d;  // d ^ d_other
    e_all = e;  // e ^ e_other

    prod = c ^ (x2 and e_all) ^ (y2 and d_all);
    if (server_num == 0)
      prod ^= (d_all and e_all);

    carry ^= prod;
  }

  return carry;
}

#endif
