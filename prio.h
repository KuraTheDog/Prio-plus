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
int twoOrder = 19;

fmpz_t Int_Modulus;
fmpz_t Int_Gen;
flint_rand_t seed;
fmpz_t *roots = nullptr, *invroots = nullptr;

typedef struct client_packet*  ClientPacket;

void SplitShare(fmpz_t val, fmpz_t A, fmpz_t B);

void init_constants(); 

void clear_constants();

#endif