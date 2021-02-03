#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <gmpxx.h>

#include <string>

extern "C" {
  #include "flint/flint.h"
  #include "flint/fmpz.h"
};


/* Goal of these constants:
Int_Gen ^ (2 ^ twoOrder) = 1 mod IntModulus.

This is for FFT.
*/

// const std::string Int_Modulus_str = "5";
// const std::string Int_Gen_str = "2";
// const int twoOrder = 2;

// const std::string Int_Modulus_str = "61";  // 97 base 16
// const std::string Int_Gen_str = "8";       // 8 base 16
// const int twoOrder = 4;

// const std::string Int_Modulus_str = "161";  // 353 base 16
// const std::string Int_Gen_str = "6";
// const int twoOrder = 5;

/* 20 ^ (2 ^ 6) = 1 mod 577 */

// const std::string Int_Modulus_str = "2681";  // 9857 base 16
// const std::string Int_Gen_str = "3a";        // 58 base 16
// const int twoOrder = 7;

// const std::string Int_Modulus_str = "10001";  // 65537 base 16
// const std::string Int_Gen_str = "2";         // p = 2^16 + 1
// const int twoOrder = 16;

// const std::string Int_Modulus_str = "8008001";  // 134250497 base 16
// const std::string Int_Gen_str = "5e1298a";      // 2^(p-1 / 2^15)
// const int twoOrder = 16;

// const std::string Int_Modulus_str = "800008001";  // 36 bit modulus
// const std::string Int_Gen_str = "10fc3989c";      // 2^(p-1 / 2^15)
// const int twoOrder = 15;

// const std::string Int_Modulus_str = "80000000080001";  // 55 bit modulus
// const std::string Int_Gen_str = "66ac804179e072";      // 2^(p-1 / 2^19)
// const int twoOrder = 19;

// const std::string Int_Modulus_str = "8000000000080001";  // 63 bit modulus
// const std::string Int_Gen_str = "22855fdf11374225";
// const int twoOrder = 19;

const std::string Int_Modulus_str = "8000000000000000080001";  // 87 bit modulus
const std::string Int_Gen_str = "2597c14f48d5b65ed8dcca";      // 17567 ^ (p-1 / 2^19)
const int twoOrder = 19;

// const std::string Int_Modulus_str = "80000000000000000000080001";  // 102 bit modulus
// const std::string Int_Gen_str = "71a9f9595f292cfd55e4c5254e";
// const int twoOrder = 19;

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

void init_roots(const int N);

#endif
