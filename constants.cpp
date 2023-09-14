#include "constants.h"

#include <iostream>

#include "fmpz_utils.h"

fmpz_t Int_Modulus;
fmpz_mod_ctx_t mod_ctx;
fmpz_t Int_Gen;
size_t mod_size;
flint_rand_t seed;
size_t nbits_mod;

void init_constants() {
  fmpz_t tmp; fmpz_init(tmp);

  fmpz_set_str(tmp, Int_Modulus_str.c_str(), 16);
  init_set_fmpz_readonly(Int_Modulus, tmp);

  fmpz_mod_ctx_init(mod_ctx, Int_Modulus);

  mod_size = fmpz_size(Int_Modulus);
  nbits_mod = fmpz_clog_ui(Int_Modulus, 2);

  fmpz_set_str(tmp, Int_Gen_str.c_str(), 16);
  init_set_fmpz_readonly(Int_Gen, tmp);

  fmpz_clear(tmp);

  std::cout << "Init constants:\n";
  std::cout << "  Int_Modulus = " << fmpz_get_ui(Int_Modulus) << "\n";
  std::cout << "    # bits: " << nbits_mod << std::endl;
  std::cout << "  Int_Gen = " << fmpz_get_ui(Int_Gen) << "\n";
  std::cout << "  twoOrder = " << twoOrder << std::endl;

  flint_randinit(seed);
}

void clear_constants() {
  flint_randclear(seed);
  fmpz_clear_readonly(Int_Modulus);
  fmpz_mod_ctx_clear(mod_ctx);
  fmpz_clear_readonly(Int_Gen);
}
