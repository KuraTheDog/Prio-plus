#include "hash.h"

#include <iostream>

#include "fmpz_utils.h"

HashStore::HashStore(
    const size_t num_hashes, const size_t input_bits, const size_t hash_range,
    flint_rand_t hash_seed_arg)
: d(num_hashes)
, l_bits(input_bits)
{
  fmpz_init_set_ui(w, hash_range);
  fmpz_init_set_ui(l, 1ULL << l_bits);
  hash_seed[0] = hash_seed_arg[0];

  // TODO: LOG2 (floor), double check if this gives floor vs ceiling.
  w_bits = 1;
  size_t w_tmp = hash_range - 1;
  while(w_tmp >>= 1)
    w_bits++;

  std::cout << "Made hash store d: " << d << ", l: 2^" << l_bits << " -> w: " << hash_range << " (" << w_bits << " bits)" << std::endl;

  // We assume w <= L. w > L striaghtforward to do, but currently not needed
  if (fmpz_cmp(w, l) > 0) {
    std::cout << "Warning: currently only supports shrinking hashes (w <= L)" << std::endl;
  }
  // assume wd < L, else no space is saved and freq is easier
  // fmpz_t wd; fmpz_init_set(wd, w); fmpz_mul_ui(wd, wd, d);
  // if (fmpz_cmp(wd, l) >= 0) {
  //   std::cout << "Warning: wd >= L, so no space saved vs standard frequency" << std::endl;
  // }
}

void HashStore::print() const {
  for (unsigned int i = 0; i < d; i++)
    print_hash(i);
}

HashStorePoly::HashStorePoly(
    const size_t num_hashes, const size_t input_bits, const size_t hash_range,
    flint_rand_t hash_seed_arg, const size_t independence)
: HashStore(num_hashes, input_bits, hash_range, hash_seed_arg)
, degree(independence - 1)
, coeff(new fmpz_t*[d])
{
  // Constant term can be mod w, rather than mod L
  // Can't make it all mod W, since e.g. (7%5)%2 = 0 != 7 % 2
  // So with multiplication involving input x, things can happen.
  // But constant term is fine, since it's a final addition of pure random.
  for (unsigned int i = 0; i < d; i++) {
    new_fmpz_array(&coeff[i], degree + 1);
    for (unsigned int j = 0; j <= degree; j++) {
      fmpz_randm(coeff[i][j], hash_seed, j == 0 ? w : l);
    }
  }
}

void HashStorePoly::eval(const unsigned int i, const unsigned int x, fmpz_t out) const {
  fmpz_zero(out);
  // f -> x(f + c_j)
  for (unsigned int j = degree; j > 0; j--) {
    fmpz_add(out, out, coeff[i][j]);
    fmpz_mul_ui(out, out, x);
    fmpz_mod(out, out, l);
  }
  // right shift, take first w_bits bits
  fmpz_fdiv_q_2exp(out, out, l_bits - w_bits);

  fmpz_add(out, out, coeff[i][0]);
  fmpz_mod(out, out, w);
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

