/* Hash family manager

Makes num_hashes hashes [input_range] -> [output_range] (parameterized by bits).

Poly:
  k-wise independent via deg k-1 polynomials
Makes a set of num_hashes random k-wise indepndent hashes via degree k-1 polynomials
Takes in a seed, so that different stores can be synced.

X-wise independent: degree x-1 polynomials

Ideally, output_range is field (p^k), especially 2^k, for perfect pairwise.
For efficiency, fine for output_range smaller and just have slightly imbalanced.
Instead treats as though 2^k for next power of 2, then mod output_range.
Does poly(x), take first k bits, then mod output_range

Could be more efficient to use fmpz_t poly, but more work.
For degree 1, efficient enough as is
*/

#ifndef PRIO_HASH_H
#define PRIO_HASH_H

#include <iostream>

#include "fmpz_utils.h"

class HashStore {
protected:
  const size_t num_hashes;
  const size_t input_bits;
  fmpz_t input_range;
  fmpz_t output_range;
  const size_t output_bits;
  flint_rand_t hash_seed;   // synced seed for consistent random hashes

public:

  HashStore(
      const size_t num_hashes, const size_t input_bits, const size_t hash_range,
      flint_rand_t hash_seed_arg);

  ~HashStore() {
    fmpz_clear(output_range);
    fmpz_clear(input_range);
    flint_randclear(hash_seed);
  }

  // Run hash i on x, returning out
  virtual void eval(const unsigned int i, const unsigned int x, fmpz_t out) const = 0;
  virtual void print_hash(const unsigned int i) const = 0;
  void print() const;
};

class HashStorePoly : public HashStore {
  const size_t degree;  // degree+1 -wise independent.

  fmpz_t** const coeff;

public:

  HashStorePoly(
      const size_t num_hashes, const size_t input_bits, const size_t hash_range,
      flint_rand_t hash_seed_arg, const size_t independence = 2);

  ~HashStorePoly() {
    for (unsigned int i = 0; i < num_hashes; i++)
      clear_fmpz_array(coeff[i], degree + 1);
    delete[] coeff;
  }

  void eval(const unsigned int i, const unsigned int x, fmpz_t out) const;
  void print_hash(const unsigned int i) const;
};

#endif