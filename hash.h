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

extern "C" {
    #include "flint/fmpz_mod_mat.h"
}

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

  virtual ~HashStore() {
    fmpz_clear(output_range);
    fmpz_clear(input_range);
    flint_randclear(hash_seed);
  }

  // Run hash i on x, returning out
  virtual void eval(const unsigned int i, const uint64_t x, fmpz_t out) const = 0;
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
    if (output_bits > 0)
      for (unsigned int i = 0; i < num_hashes; i++)
        clear_fmpz_array(coeff[i], degree + 1);
    delete[] coeff;
  }

  void eval(const unsigned int i, const uint64_t x, fmpz_t out) const;
  void print_hash(const unsigned int i) const;
};

class HashStoreBit : public HashStore {
  // How big groups of hashes are that are all on the same value, for solving.
  const size_t group_size;

  // General "input dimension" = input_bits + 1
  // Constant since otherwise h(0) = 0 always.
  const size_t dim;
  fmpz_mod_mat_t coeff;

public:

  // Default group size 0 -> num_hashes, aka 1 group.
  HashStoreBit(
      const size_t num_hashes, const size_t input_bits, const size_t hash_range,
      flint_rand_t hash_seed_arg, const size_t group_size_arg = 0);

  ~HashStoreBit() {
    fmpz_mod_mat_clear(coeff);
  };

  void eval(const unsigned int i, const uint64_t x, fmpz_t out) const;
  void print_hash(const unsigned int i) const;
  void print_coeff() const {fmpz_mod_mat_print_pretty(coeff);}

  // |values| = group_size
  // Uses Group(group_num) of hashes.
  // Uses first [dim] of the hashes and values to solve AX=B
  // Uses extras [dim -> group size] for validation.
  //   Returns number of INVALID checks. 0 is best
  int solve(const unsigned int group_num, const fmpz_t* const values, uint64_t& ans) const;
};

#endif
