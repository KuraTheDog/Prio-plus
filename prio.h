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

const std::string Int_Modulus_str = "8000000000000000080001";
const std::string Int_Gen_str = "2597c14f48d5b65ed8dcca";
const int twoOrder = 19;

extern fmpz_t Int_Modulus;
extern fmpz_t Int_Gen;
extern flint_rand_t seed;
extern fmpz_t *roots, *invroots;

void init_constants(); 

void clear_constants();

#endif