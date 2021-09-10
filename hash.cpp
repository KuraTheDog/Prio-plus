/* Only need pairwise independnet hashes for single countmin

X-wise independent: degree x-1 polynomials
We just need 2 for count-min. 
TODO: double check for heavy that only need 2

Could be more efficient to use fmpz_t poly, but more work. 
For degree 1, easy enough to use 1

*/

#include "hash.h"

#include <iostream>

#include "fmpz_utils.h"

void HashStore::eval(const unsigned int i, const fmpz_t x, fmpz_t out) {
  // Currently hard code for degree 1
  fmpz_set(out, coeff[i][0]);
  fmpz_addmul(out, coeff[i][1], x);
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
