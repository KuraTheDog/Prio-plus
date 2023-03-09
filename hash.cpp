#include "hash.h"

#include <iostream>

#include "fmpz_utils.h"
#include "utils.h"

HashStore::HashStore(
    const size_t num_hashes, const size_t input_bits, const size_t hash_range,
    flint_rand_t hash_seed_arg)
: num_hashes(num_hashes)
, input_bits(input_bits)
, output_bits(hash_range > 1 ? LOG2(hash_range - 1) : 0)
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

HashStoreBit::HashStoreBit(
    const size_t num_hashes, const size_t input_bits, const size_t hash_range,
    flint_rand_t hash_seed_arg, const size_t group_size_arg)
: HashStore(num_hashes, input_bits, hash_range, hash_seed_arg)
, group_size(group_size_arg ? group_size_arg : num_hashes)
, dim(input_bits + 1)
{
  const size_t num_groups = num_hashes / group_size;

  // std::cout << "  Bit store: " << num_groups << " groups of size " << group_size << ", ";
  // std::cout << "  dim : " << dim << ", validate: " << group_size - dim << std::endl;

  if (num_hashes % group_size != 0) {
    std::cout << "Warning: Can't evenly split " << num_hashes << " hashes into size " << group_size << " groups." << std::endl;
  }
  if (group_size < dim) {
    std::cout << "Warning: Group size " << group_size << " smaller than input dim " << dim << " groups." << std::endl;
  }

  fmpz_mod_mat_init(coeff, num_hashes, dim, output_range);
  // Random for all/extra outside of main solvable square
  fmpz_mod_mat_randtest(coeff, hash_seed);

  // Goal: For each group, first square bit is singular, for inversion purposes
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
  fmpz_add(out, out, fmpz_mod_mat_entry(coeff, i, input_bits));
}

void HashStoreBit::print_hash(const unsigned int i) const {
  std::cout << "Hash " << i << ": ";
  for (unsigned int j = 0; j < dim; j++) {
    fmpz_print(fmpz_mod_mat_entry(coeff, i, j));
    if (j < input_bits) std::cout << " * x_" << j << " + ";
  }
  std::cout << std::endl;
}

int HashStoreBit::solve(const unsigned int group_num, const fmpz_t* const values, uint64_t& ans) const {
  const size_t start = group_num * group_size;

  // Init
  fmpz_mod_mat_t A;
  fmpz_mod_mat_t B; fmpz_mod_mat_init(B, dim, 1, output_range);
  fmpz_mod_mat_t X; fmpz_mod_mat_init(X, dim, 1, output_range);
  for (unsigned int i = 0; i < dim; i++) {
    fmpz_set(fmpz_mod_mat_entry(B, i, 0), values[i]);
  }
  fmpz_mod_mat_window_init(A, coeff, start, 0, start + dim, dim);

  // Solve
  int ret = fmpz_mod_mat_solve(X, A, B);
  fmpz_mod_mat_window_clear(A);
  fmpz_mod_mat_clear(B);
  if (ret == 0) {
    // std::cout << "WARNING: Somehow couldn't solve." << std::endl;
    fmpz_mod_mat_clear(X);
    // Case should not happen
    return 100001;
  }

  // Extract
  ans = 0;
  for (unsigned int i = 0; i < input_bits; i++) {
    if (fmpz_is_one(fmpz_mod_mat_entry(X, i, 0))) {
      // fmpz_add_ui(ans, ans, 1ULL << i);
      ans += (1ULL << i);
    } else if (fmpz_is_zero(fmpz_mod_mat_entry(X, i, 0))) {
      continue;
    } else {
      // std::cout << "WARNING: Solution has non-bit value" << std::endl;
      fmpz_mod_mat_clear(X);
      // Case likely because inconsistent in the first n+1.
      return 100002;
    }
  }
  if (not fmpz_is_one(fmpz_mod_mat_entry(X, input_bits, 0))) {
    // std::cout << "WARNING: Solution has non-1 constant term" << std::endl;
    fmpz_mod_mat_clear(X);
    // Inconsistent in first n+1, or pure empty.
    return 100003;
  }
  fmpz_mod_mat_clear(X);

  // Check
  // TODO: it would be nicer if it could do least error version:
  //   E.g. if 2n equations when need n, and 1 bad, it can solve on the 2n-1 good.
  //      (vs including bad makes many more invalid)
  //   Solve / inverse don't seem to like non-square matricies
  //   and least squares seems nontrivial over modular context

  int num_invalid = 0;
  fmpz_t tmp; fmpz_init(tmp);
  for (unsigned int i = input_bits; i < group_size; i++) {
    // Ensure hash_{start + i}(ans) = values[i]
    eval(start + i, ans, tmp);
    if (not fmpz_equal(tmp, values[i])) {
      num_invalid++;
    }
  }
  fmpz_clear(tmp);

  return num_invalid;
}
