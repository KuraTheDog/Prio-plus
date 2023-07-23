#include "hash.h"

#include <iostream>

#include "fmpz_utils.h"
#include "utils.h"

HashStore::HashStore(
    const size_t num_hashes, const size_t input_bits, const size_t hash_range,
    flint_rand_t hash_seed_arg)
: num_hashes(num_hashes)
, input_bits(input_bits)
, output_bits(LOG2(hash_range))
{
  fmpz_init_set_ui(output_range, hash_range);
  fmpz_init_set_ui(input_range, 1ULL << input_bits);
  hash_seed[0] = hash_seed_arg[0];

  // std::cout << "Made hash store num_hashes: " << num_hashes << ", input_range: 2^" << input_bits;
  // std::cout << " -> output_range: " << hash_range << " (" << output_bits << " bits)" << std::endl;
}

void HashStore::print() const {
  for (unsigned int i = 0; i < num_hashes; i++)
    print_hash(i);
}

HashStorePoly::HashStorePoly(
    const size_t num_hashes, const size_t input_bits, const size_t hash_range,
    flint_rand_t hash_seed_arg, const size_t independence)
: HashStore(num_hashes, input_bits, hash_range, hash_seed_arg)
, degree(independence - 1)
, coeff(new fmpz_t*[num_hashes])
{
  // Input [input_range], output [output_range]. What modulo can coefficients be?
  // Final answer is always modulo output_range. Can coefficients be pre-mod output_range?
  // e.g. (7%5)%2 = 0 != 7 % 2
  // Constants fixed, but what about when multiplied by constant?
  // Want distribution same.
  // So with multiplication involving input x, things can happen.
  // But constant term should be fine, since it's a final addition of pure random.
  // We also just want "good enough" distribution.
  // Also, it's not quite polynomial since shift magic.
  // std::cout << "  Poly store degree: " << degree << std::endl;
  // Can also run into coprime issues. E.g. 3 bits in and out, 0 + 4x, always outputs 0 or 4
  if (output_bits == 0) return;
  for (unsigned int i = 0; i < num_hashes; i++) {
    new_fmpz_array(&coeff[i], degree + 1);
    // Enforce not constant. Only really issue for small input ranges, and low degree
    bool run = true;
    while(run) {
      for (unsigned int j = 0; j <= degree; j++) {
        fmpz_randm(coeff[i][j], hash_seed, j == 0 ? output_range : input_range);
      }
      for (unsigned int j = 1; j <= degree; j++) {
        if (not fmpz_is_zero(coeff[i][j])) {
          run = false;
          break;
        }
      }
    }
  }
}

void HashStorePoly::eval(const unsigned int i, const uint64_t x, fmpz_t out) const {
  fmpz_zero(out);
  if (output_bits == 0) return;
  // f -> x(f + c_j)
  for (unsigned int j = degree; j > 0; j--) {
    fmpz_add(out, out, coeff[i][j]);
    fmpz_mul_ui(out, out, x);
    fmpz_mod(out, out, input_range);
  }
  // right shift, take first output_bits bits
  if (input_bits > output_bits)
    fmpz_fdiv_q_2exp(out, out, input_bits - output_bits);

  fmpz_add(out, out, coeff[i][0]);
  fmpz_mod(out, out, output_range);
}

void HashStorePoly::print_hash(const unsigned int i) const {
  std::cout << "Hash " << i << ": ";
  if (output_bits == 0) {
    std::cout << "0" << std::endl;
    return;
  }
  for (unsigned int j = 0; j <= degree; j++) {
    fmpz_print(coeff[i][j]);
    if (j > 0) std::cout << " x";
    if (j > 1) std::cout << "^" << j;
    if (j < degree) std::cout << " + ";
  }
  std::cout << std::endl;
}

void HashStoreBit::build_coefficients() {
  fmpz_mod_mat_init(coeff, num_hashes, dim, output_range);

  // Random for all/extra outside of main solvable square
  fmpz_mod_mat_randtest(coeff, hash_seed);

  // Goal: For each group, first square bit is singular, for inversion purposes
  inverses = new fmpz_mod_mat_t[num_groups];
  fmpz_mod_mat_t window;
  for (unsigned int i = 0; i < num_groups; i++) {
    const size_t start = i * group_size;
    // By reference element access
    fmpz_mod_mat_window_init(window, coeff, start, 0, start + dim, dim);

    // Minimal setting, so need to randomize after
    fmpz_mod_mat_randrank(window, hash_seed, dim);
    // How much to get dense random? Flint t-det.c uses 2n+1, t-solve.c uses n*n+1
    // Also Flint tests prefer random: 1+n_randint(hash_seed, 2 * input_bits + 1)
    // a) which seems like it'd bias more towards more 0's?
    // b) Docs already say "at most count", so seems already random amount.
    fmpz_mod_mat_randops(window, 4 * dim + 1, hash_seed);

    // Also compute inverse ahead of time
    // By above, always convertable
    if (is_solving) {
      fmpz_mod_mat_init(inverses[i], dim, dim, output_range);
      fmpz_mod_mat_inv(inverses[i], window);
      // _fmpz_mod_mat_set_mod(inverses[i], Int_Modulus);
    }
    fmpz_mod_mat_window_clear(window);
  }
}

void HashStoreBit::eval(const unsigned int i, const uint64_t x, fmpz_t out) const {
  fmpz_zero(out);
  for (unsigned int j = 0; j < input_bits; j++) {
    // if (fmpz_tstbit(x, j)) {
    if ( (x>>j)%2 == 1 ) {
      fmpz_add(out, out, fmpz_mod_mat_entry(coeff, i, j));
      fmpz_mod(out, out, output_range);
    }
  }
  // Constant
  if (inconsistency_solving)
    fmpz_add(out, out, fmpz_mod_mat_entry(coeff, i, input_bits));
}

void HashStoreBit::print_hash(const unsigned int i) const {
  std::cout << "Hash " << i << ": ";
  for (unsigned int j = 0; j < dim; j++) {
    fmpz_print(fmpz_mod_mat_entry(coeff, i, j));
    if (j < input_bits) {
      std::cout << " * x_" << j;
      if (inconsistency_solving or j < (input_bits - 1))
        std::cout << " + ";
    }
  }
  std::cout << std::endl;
}

void HashStoreBit::solve_shares(const unsigned int group_num,
    const fmpz_t* const values, fmpz_t* const ans) const {
  if (!is_solving) {
    error_exit("Hash store trying to solve when is_solving is false");
  }

  // Coeff * X = values -> X = coeff^-1 * values
  // TODO: better modulus, e.g. if values is shared (int modulus)
  fmpz_mod_mat_t V; fmpz_mod_mat_init(V, dim, 1, output_range);
  // TODO: X's range should be input? But mod alignment weird.
  fmpz_mod_mat_t X; fmpz_mod_mat_init(X, dim, 1, output_range);

  for (unsigned int i = 0; i < dim; i++) {
    fmpz_set(fmpz_mod_mat_entry(V, i, 0), values[i]);
  }

  fmpz_mod_mat_mul(X, inverses[group_num], V);

  fmpz_mod_mat_clear(V);

  for (unsigned int i = 0; i < dim; i++) {
    fmpz_set(ans[i], fmpz_mod_mat_entry(X, i, 0));
    // std::cout << "X[" << i << "] = " << fmpz_get_ui(ans[i]) << std::endl;
  }
  fmpz_mod_mat_clear(X);
}

int HashStoreBit::solve_extract(const unsigned int group_num, const fmpz_t* const values,
    const fmpz_t* const vec, uint64_t& ans) const {
  ans = 0;
  // Extract answer
  for (unsigned int i = 0; i < input_bits; i++) {
    if (fmpz_is_one(vec[i])) {
      // fmpz_add_ui(ans, ans, 1ULL << i);
      ans += (1ULL << i);
    } else if (fmpz_is_zero(vec[i])) {
      continue;
    } else {
      // std::cout << "WARNING: Solution has non-bit value" << std::endl;
      // Case occurs when inconsistent in the first n+1, which is not unlikely
      return 100002;
    }
  }

  // Check constant term
  if (inconsistency_solving and
      not fmpz_is_one(vec[input_bits])) {
    // std::cout << "WARNING: Solution has non-1 constant term" << std::endl;
    // Inconsistent in first n+1, or pure empty.
    return 100003;
  }

  // Compare against remaing values
  const size_t offset = group_num * group_size;
  int num_invalid = 0;
  if (inconsistency_solving) {
    fmpz_t tmp; fmpz_init(tmp);
    for (unsigned int i = input_bits; i < group_size; i++) {
      // Ensure hash_{offset + i}(ans) = values[i]
      eval(offset + i, ans, tmp);
      if (not fmpz_equal(tmp, values[i])) {
        num_invalid++;
      }
    }
    fmpz_clear(tmp);
  }

  return num_invalid;
}

int HashStoreBit::solve(const unsigned int group_num,
    const fmpz_t* const values, uint64_t& ans) const {
  fmpz_t* X; new_fmpz_array(&X, dim);
  solve_shares(group_num, values, X);

  // Extract
  int ret = solve_extract(group_num, values, X, ans);
  clear_fmpz_array(X, dim);

  return ret;
}
