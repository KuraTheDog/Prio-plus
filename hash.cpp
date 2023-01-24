#include "hash.h"

#include <iostream>

#include "fmpz_utils.h"
#include "utils.h"

HashStore::HashStore(
    const size_t num_hashes, const size_t input_bits, const size_t hash_range,
    flint_rand_t hash_seed_arg)
: num_hashes(num_hashes)
, input_bits(input_bits)
, output_bits(LOG2(hash_range - 1))
{
  fmpz_init_set_ui(output_range, hash_range);
  fmpz_init_set_ui(input_range, 1ULL << input_bits);
  hash_seed[0] = hash_seed_arg[0];

  std::cout << "Made hash store num_hashes: " << num_hashes << ", input_range: 2^" << input_bits << " -> output_range: " << hash_range << " (" << output_bits << " bits)" << std::endl;

  // We assume output_range <= input_range. output_range > input_range striaghtforward to do, but currently not needed
  if (fmpz_cmp(output_range, input_range) > 0) {
    std::cout << "Warning: currently only supports shrinking hashes (output_range <= input_range)" << std::endl;
  }
  // assume wd < input_range, else no space is saved and freq is easier
  // fmpz_t wd; fmpz_init_set(wd, output_range); fmpz_mul_ui(wd, wd, num_hashes);
  // if (fmpz_cmp(wd, input_range) >= 0) {
  //   std::cout << "Warning: wd >= input_range, so no space saved vs standard frequency" << std::endl;
  // }
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
  for (unsigned int i = 0; i < num_hashes; i++) {
    new_fmpz_array(&coeff[i], degree + 1);
    for (unsigned int j = 0; j <= degree; j++) {
      fmpz_randm(coeff[i][j], hash_seed, j == 0 ? output_range : input_range);
    }
  }
}

void HashStorePoly::eval(const unsigned int i, const unsigned int x, fmpz_t out) const {
  fmpz_zero(out);
  // f -> x(f + c_j)
  for (unsigned int j = degree; j > 0; j--) {
    fmpz_add(out, out, coeff[i][j]);
    fmpz_mul_ui(out, out, x);
    fmpz_mod(out, out, input_range);
  }
  // right shift, take first output_bits bits
  fmpz_fdiv_q_2exp(out, out, input_bits - output_bits);

  fmpz_add(out, out, coeff[i][0]);
  fmpz_mod(out, out, output_range);
}

void HashStorePoly::print_hash(const unsigned int i) const {
  std::cout << "Hash " << i << ": ";
  for (unsigned int j = 0; j <= degree; j++) {
    fmpz_print(coeff[i][j]);
    if (j > 0) std::cout << " x";
    if (j > 1) std::cout << "^" << j;
    if (j < degree) std::cout << " + ";
  }
  std::cout << std::endl;
}

