#include "constants.h"

#include <iostream>

#include "fmpz_utils.h"

fmpz_t Int_Modulus;
fmpz_t Int_Gen;
flint_rand_t seed;

void init_constants() {
    fmpz_init(Int_Modulus);
    fmpz_set_str(Int_Modulus,Int_Modulus_str.c_str(),16);
    fmpz_init(Int_Gen);
    fmpz_set_str(Int_Gen,Int_Gen_str.c_str(),16);

    std::cout << "Init constants: " << std::endl;
    std::cout << "  Int_Modulus = "; fmpz_print(Int_Modulus); std::cout << std::endl;
    std::cout << "  Int_Gen = "; fmpz_print(Int_Gen); std::cout << std::endl;
    std::cout << "  twoOrder = " << twoOrder << std::endl;
    flint_randinit(seed);
}

void clear_constants() {
    flint_randclear(seed);
    fmpz_clear(Int_Modulus);
    fmpz_clear(Int_Gen);
}
