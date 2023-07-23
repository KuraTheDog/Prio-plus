/* Hash family manager

Makes num_hashes hashes [input_range] -> [output_range] (parameterized by bits).
With stacking/resuse, num_hashes can include many groups of hashes.

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

  size_t get_input_bits() const { return input_bits; };
  size_t get_num_hashes() const { return num_hashes; };
};

class HashStorePoly : public HashStore {
protected:
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

  uint64_t get_coeff(const size_t i, const size_t j) const {
    return fmpz_get_ui(coeff[i][j]);
  }
};

// R = num_bits hashes
// survives(r, x) with probability 1 / 2^r
// "0" is always survive, up to "R" is one survivor
// Extra to force more diversity: ax+b mod 2^R, with a odd.
class HashStoreHalf : public HashStorePoly{
protected:
  const size_t R;
public:
  HashStoreHalf(const size_t num_bits, flint_rand_t hash_seed_arg)
  : HashStorePoly(num_bits, num_bits, 1ULL << num_bits, hash_seed_arg, 2)
  , R(num_bits)
  {
    // resample, since shennanigans with co-primeness. so force a = odd
    fmpz_t half_range;
    fmpz_init_set_ui(half_range, 1ULL << (R - 1));
    for (unsigned int i = 0; i < num_bits; i++) {
      fmpz_randm(coeff[i][1], hash_seed, half_range);
      fmpz_mul_ui(coeff[i][1], coeff[i][1], 2);
      fmpz_add_ui(coeff[i][1], coeff[i][1], 1);
    }
    fmpz_clear(half_range);
  }

  bool survives(const size_t i, const uint64_t x) const {
    if (i == 0)
      return true;

    fmpz_t tmp; fmpz_init(tmp);
    eval(i - 1, x, tmp);
    fmpz_fdiv_q_2exp(tmp, tmp, R - i);
    bool ret = (bool) fmpz_is_zero(tmp);
    fmpz_clear(tmp);
    return ret;
  }
};

/*
  When group_size = input_bits (default), just solves, with no inconsitency checks
  When group_size > input_bits, extra equations are checked against solved value
    Returns # of inconsistent, so 0 is correct.
    Also includes "error" returns of 10000, 10001, and 10002
    for other inconsitencies
*/
class HashStoreBit : public HashStore {
protected:
  // How big groups of hashes are that are all on the same value, for solving.
  // group size * num_groups gives num_hashes
  // group size should be >= input bits
  // If and only if it's equal, then no inconsistnecy solving
  const size_t num_groups;
  const size_t group_size;

  // If false, doesn't precompute inverses, faster but then can't solve.
  // Usefor for clients just evaluating on hash
  const bool is_solving;
  // If true, extra equations checking against solved value
  const bool inconsistency_solving;

  // General "input dimension" = input_bits + 1
  // Constant since otherwise h(0) = 0 always.
  const size_t dim;
  fmpz_mod_mat_t coeff;

  fmpz_mod_mat_t* inverses = nullptr;

public:

  HashStoreBit(
    const size_t num_groups, const size_t group_size, const size_t input_bits,
    const size_t hash_range, flint_rand_t hash_seed_arg, const bool is_solving = true)
  : HashStore(num_groups * group_size, input_bits, hash_range, hash_seed_arg)
  , num_groups(num_groups)
  , group_size(group_size)
  , is_solving(is_solving)
  , inconsistency_solving(group_size > input_bits)
  , dim(input_bits + ((int) inconsistency_solving))
  {
    // std::cout << (inconsistency_solving ? "" : "not ") << "inconsistency_solving" << std::endl;
    // std::cout << "  Bit store: " << num_groups << " groups of size " << group_size << ", ";
    // std::cout << "  dim : " << dim << ", validate: " << group_size - dim << std::endl;

    if (group_size < input_bits) {
      std::cout << "Warning: Group size " << group_size << " smaller than input bits " << input_bits << std::endl;
    }

    build_coefficients();
  }

  ~HashStoreBit() {
    fmpz_mod_mat_clear(coeff);
    if (is_solving) {
      for (unsigned int i = 0; i < num_groups; i++) {
        fmpz_mod_mat_clear(inverses[i]);
      }
      delete[] inverses;
    }
  };

  void build_coefficients();
  void eval(const unsigned int i, const uint64_t x, fmpz_t out) const;
  void print_hash(const unsigned int i) const;
  void print_coeff() const { fmpz_mod_mat_print_pretty(coeff); }
  void print_inv(const size_t i) const { fmpz_mod_mat_print_pretty(inverses[i]); }
  size_t get_num_groups() const { return num_groups; };
  uint64_t get_inv_coeff(const size_t idx, const size_t i, const size_t j) const {
    return fmpz_get_ui(fmpz_mod_mat_entry(inverses[idx], i, j));
  }

  // Call before doing solve on shared values, to shift mod context
  void adjust_mod(const fmpz_t new_mod) {
    for (unsigned int i = 0; i < num_groups; i++) {
      _fmpz_mod_mat_set_mod(inverses[i], new_mod);
    }
  }

  // |values| = group_size
  // Uses Group(group_num) of hashes.
  // Uses inverse first [dim] of the hashes and values to solve X=A^-1 B
  // Uses extras [dim -> group size] for validation.
  //   Returns number of INVALID checks. 0 is best
  int solve(const unsigned int group_num, const fmpz_t* const values, uint64_t& ans) const;
  // Finds ans, in the corresponding vector form.
  void solve_shares(const unsigned int group_num, const fmpz_t* const values, fmpz_t* const ans) const;
  // Given "vector" representation vec, extracts answer into ans
  // Also uses extras from values to validate.
  // Note that if not inconsistency solving, then values is unused
  int solve_extract(const unsigned int group_num, const fmpz_t* const values,
      const fmpz_t* const vec, uint64_t& ans) const;
};

#endif
