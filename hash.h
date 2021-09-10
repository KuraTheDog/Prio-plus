/* Hash family manager

Makes a set of d random k-wise indpeendent hashes via degree k-1 polynomials
Takes in a seed, so that different stores can be synced.
Hashes are range (mod) w, which should be a field size (p^n)

Since clients send boolean shares mod w, maximizing w by rounding up to nearest power of 2 gives better error at no additional cost.
*/

#ifndef PRIO_HASH_H
#define PRIO_HASH_H

#include <iostream>

#include "fmpz_utils.h"

class HashStore {
  fmpz_t w;  // hash range. Should be p^n. power of 2 gives best efficiency?
  const size_t d;   // # hashes
  const size_t degree = 1;  // degree+1 -wise independent. 
  flint_rand_t hash_seed;  // synced seed for same random hashes

  fmpz_t** coeff;

public: 
  HashStore(const size_t w_arg, const size_t d, flint_rand_t hash_seed_arg)
  : d(d) {
    fmpz_init(w);
    fmpz_set_ui(w, w_arg);
    hash_seed[0] = hash_seed_arg[0];
    // flint_randinit(hash_seed);

    coeff = new fmpz_t*[d];

    // are we ok with constant poly? aka a + 0x?
    // or even linear coeff makes it miss half
    for (unsigned int i = 0; i < d; i++) {
      new_fmpz_array(&coeff[i], degree + 1);
      for (unsigned int j = 0; j <= degree; j++)
        fmpz_randm(coeff[i][j], hash_seed, w);
    }
  }

  ~HashStore() {
    fmpz_clear(w);
    flint_randclear(hash_seed);
    for (unsigned int i = 0; i < d; i++)
      clear_fmpz_array(coeff[i], degree + 1);
    delete[] coeff;
  }

  // Run hash i on x, returning out
  void eval(const unsigned int i, const fmpz_t x, fmpz_t out);
  void print_hash(const unsigned int i);
};

#endif