#include "hash.h"

#include <iostream>

#include "fmpz_utils.h"

void HashStore::eval(const unsigned int i, const unsigned int x, fmpz_t out) {
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

void HashStore::print_hash(const unsigned int i) {
  std::cout << "Hash " << i << ": ";
  for (unsigned int j = 0; j <= degree; j++) {
    fmpz_print(coeff[i][j]);
    if (j > 0) std::cout << " x";
    if (j > 1) std::cout << "^" << j;
    if (j < degree) std::cout << " + ";
  }
  std::cout << std::endl;
}
