#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <gmpxx.h>

#include <string>

extern "C" {
  #include "flint/flint.h"
  #include "flint/fmpz.h"
  #include "flint/fmpz_mod.h"
};

/*
Int_Gen ^ (2 ^ twoOrder) = 1 mod IntModulus.
This is for FFT.

Int_Modulus should have p = 1 mod 2^k for some k.
This is for BFV.

Note:
g = a^((p-1) / 2^k) satisfies g^(2^k) = 1 mod p.
Smaller k' may work, so differnet a may be required.
*/

// const std::string Int_Modulus_str = "5";
// const std::string Int_Gen_str = "2";
// const int twoOrder = 2;

// const std::string Int_Modulus_str = "61";  // 97 base 16
// const std::string Int_Gen_str = "8";       // 8 base 16
// const int twoOrder = 4;

// const std::string Int_Modulus_str = "101";  // 257 base 16
// const std::string Int_Gen_str = "3";
// const int twoOrder = 8;

// const std::string Int_Modulus_str = "3401";  // 13313, 14 bit modulus
// const std::string Int_Gen_str = "3";         // 3^(p-1 / 2^10)
// const int twoOrder = 10;

// const std::string Int_Modulus_str = "10001";  // 65537, 17 bits modulus
// const std::string Int_Gen_str = "3";          // p = 2^16 + 1
// const int twoOrder = 16;

// const std::string Int_Modulus_str = "8008001";  // 28 bit modulus
// const std::string Int_Gen_str = "1CF4C77";      // 3^(p-1 / 2^15)
// const int twoOrder = 15;

// const std::string Int_Modulus_str = "800008001";  // 36 bit modulus
// const std::string Int_Gen_str = "7946479F9";      // 3^(p-1 / 2^15)
// const int twoOrder = 15;

const std::string Int_Modulus_str = "800006880001";  // 48 bit modulus
const std::string Int_Gen_str = "3503101C8855";        // 7^(p-1 / 2^19)
const int twoOrder = 19;

// Below here doesn't work with PALISADE triples

// const std::string Int_Modulus_str = "80000000080001";  // 55 bit modulus
// const std::string Int_Gen_str = "4359077260C2D6";      // 3^(p-1 / 2^19)
// const int twoOrder = 19;

// const std::string Int_Modulus_str = "8000000000080001";  // 63 bit modulus
// const std::string Int_Gen_str = "22855fdf11374225";      // 965081^(p-1 / 2^19)
// const int twoOrder = 19;

// const std::string Int_Modulus_str = "8000000000000000080001";  // 87 bit modulus
// const std::string Int_Gen_str = "2597c14f48d5b65ed8dcca";      // 17567 ^ (p-1 / 2^19)
// const int twoOrder = 19;

// const std::string Int_Modulus_str = "80000000000000000000080001";  // 102 bit modulus
// const std::string Int_Gen_str = "71a9f9595f292cfd55e4c5254e";
// const int twoOrder = 19;

extern fmpz_t Int_Modulus;        // Large prime modulus
extern fmpz_mod_ctx_t mod_ctx;    // Modulus context
extern fmpz_t Int_Gen;            // Generates subgroup order 2^twoOrder in Zp
extern flint_rand_t seed;         // Global random seed, for fmpz_randm, etc.
extern size_t mod_size;           // for net share, fmpz_size(Int_Modulus);

/* Nth roots of unity, used for FFT.
   over 2^twoOrder.
   N is a power of 2, 2^k.
   g is Int_Gen
   r is an Nth root of unity. Done by 2^(ord - k).
Then we have roots = 1, r, r^2, ..., r^{N-1}, and their inverses.
root2 are the 2Nth roots of unity.
*/
extern size_t num_roots;  // Tracker for N, so we can clear and rebuild
extern fmpz_t *roots, *invroots, *roots2;

void init_constants();

void clear_constants();

void init_roots(const size_t N);

#endif
