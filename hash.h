/* Hash family manager

Makes a set of d random k-wise indpeendent hashes via degree k-1 polynomials
Takes in a seed, so that different stores can be synced.

X-wise independent: degree x-1 polynomials
We just need 2 (pairwise) for count-min.

Ideally, w is field (p^k), especially 2^k, for perfect pairwise.
For efficiency, fine for w smaller and just have slightly imbalanced. Instead treats as though 2^k for next power of 2, then mod w.
Does poly(x), take first k bits, then mod w

Could be more efficient to use fmpz_t poly, but more work. 
For degree 1, efficient enough as is
*/

#ifndef PRIO_HASH_H
#define PRIO_HASH_H

#include <iostream>

#include "fmpz_utils.h"

class HashStore {
  const size_t d;   // # hashes
  const size_t l_bits;  // input bits
  fmpz_t l;  // input range
  fmpz_t w;  // hash range
  size_t w_bits;  // # bits to store w
  const size_t degree;  // degree+1 -wise independent. 
  flint_rand_t hash_seed;  // synced seed for same random hashes

  fmpz_t** coeff;

public: 
  HashStore(const size_t d, const size_t l_bits, const size_t w_arg,
            flint_rand_t hash_seed_arg, const size_t independence = 2)
  : d(d)
  , l_bits(l_bits)
  , degree(independence - 1)
  {
    fmpz_init(w); fmpz_set_ui(w, w_arg);
    fmpz_init(l); fmpz_set_ui(l, 1ULL << l_bits);
    hash_seed[0] = hash_seed_arg[0];

    w_bits = 1;
    size_t w_tmp = w_arg - 1;
    while(w_tmp >>= 1)
      w_bits++;

    // std::cout << "Hash store d: " << d << ", l: 2^" << l_bits << " -> w: " << w_arg << " (" << w_bits << " bits)" << std::endl;

    // TODO: sanity checks:
    // We assume w <= L. w > L striaghtforward to do, but currently not needed
    if (fmpz_cmp(w, l) > 0) {
      std::cout << "Warning: currently only supports shrinking hashes (w <= L)" << std::endl;
    }
    // assume wd < L, else no space is saved and freq is easier
    fmpz_t wd; fmpz_init_set(wd, w); fmpz_mul_ui(wd, wd, d);
    if (fmpz_cmp(wd, l) >= 0) {
      std::cout << "Warning: wd >= L, so no space saved vs standard frequency" << std::endl;
    }

    coeff = new fmpz_t*[d];
    // Constant term can be mod w, rather than mod L
    for (unsigned int i = 0; i < d; i++) {
      new_fmpz_array(&coeff[i], degree + 1);
      for (unsigned int j = 0; j <= degree; j++) {
        if (j == 0) {
          fmpz_randm(coeff[i][j], hash_seed, w);
        } else {
          fmpz_randm(coeff[i][j], hash_seed, l);
        }
      }
    }
  }

  ~HashStore() {
    fmpz_clear(w);
    fmpz_clear(l);
    flint_randclear(hash_seed);
    for (unsigned int i = 0; i < d; i++)
      clear_fmpz_array(coeff[i], degree + 1);
    delete[] coeff;
  }

  // Run hash i on x, returning out
  void eval(const unsigned int i, const unsigned int x, fmpz_t out) const;
  void print_hash(const unsigned int i) const;
};

#endif