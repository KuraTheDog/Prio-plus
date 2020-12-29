#ifndef PRIO_H
#define PRIO_H


#include <vector>
#include <string>
#include <iostream>
#include <gmpxx.h>

extern "C" {
  #include "flint/flint.h"
  #include "flint/fmpz.h"
};

/* Simpler testing constants.
   58 ^ (2 ^ 7) = 1 mod 9857 */
// const std::string Int_Modulus_str = "2681";  // 9857 base 16
// const std::string Int_Gen_str = "3a";  // 58 base 16
// const int twoOrder = 7;

// const std::string Int_Modulus_str = "8000000000080001";  // 63 bit modulus
// const std::string Int_Gen_str = "22855fdf11374225";
// const int twoOrder = 19;          // We have group order 2^twoOrder, for use in FFT.

const std::string Int_Modulus_str = "8000000000000000080001";  // 87 bit modulus
const std::string Int_Gen_str = "2597c14f48d5b65ed8dcca";
const int twoOrder = 19;          // We have group order 2^twoOrder, for use in FFT.

// const std::string Int_Modulus_str = "80000000000000000000080001";  // 102 bit modulus
// const std::string Int_Gen_str = "71a9f9595f292cfd55e4c5254e";
// const int twoOrder = 19;          // We have group order 2^twoOrder, for use in FFT.

extern fmpz_t Int_Modulus;        // Large prime modulus
extern fmpz_t Int_Gen;            // Generates subgroup order 2^twoOrder in Zp
extern flint_rand_t seed;         // Global random seed, for fmpz_randm, etc.

/* Nth roots of unity, used for FFT.
   over 2^twoOrder.
   N is a power of 2, 2^k.
   g is Int_Gen
   r is an Nth root of unity. Done by 2^(ord - k).
Then we have roots = 1, r, r^2, ..., r^{N-1}, and their inverses.
root2 are the 2Nth roots of unity.
TODO: This is for fixed N, so caching this would be better.
*/
extern fmpz_t *roots, *invroots, *roots2;

void init_constants(); 

void clear_constants();

#endif
