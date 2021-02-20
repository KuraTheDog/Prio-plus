#include "constants.h"

#include <iostream>

#include "fmpz_utils.h"

fmpz_t Int_Modulus;
fmpz_t Int_Gen;
flint_rand_t seed;
fmpz_t *roots = nullptr, *invroots = nullptr, *roots2 = nullptr;
size_t num_roots;

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

void init_roots(const int N) {
    if (num_roots == N) {
        return;
    }
    if (roots != nullptr) {
        clear_fmpz_array(roots, num_roots);
        clear_fmpz_array(invroots, num_roots);
        clear_fmpz_array(roots2, 2 * num_roots);
    }

    num_roots = N;

    new_fmpz_array(&roots, N);
    new_fmpz_array(&invroots, N);
    new_fmpz_array(&roots2, 2 * N);

    int step_size = (1 << twoOrder)/N;  // 2^(twoOrder - log_2 N)
    fmpz_t g_, ginv_, ghalf_;
    fmpz_init(g_);      // Generator (Int_Gen)
    fmpz_init(ginv_);   // Inverse of g_
    fmpz_init(ghalf_);  // g_^(step/2), for 2N roots.

    fmpz_invmod(ginv_,Int_Gen,Int_Modulus);

    /*
    N = 2^k, so stepsize = 2^(Ord - k).
    g_ = gen^stepsize.
    So g_^N = gen^(2^ord) = 1, by fermat little.
    */
    fmpz_powm_ui(g_,Int_Gen,step_size, Int_Modulus);
    fmpz_powm_ui(ginv_,ginv_,step_size, Int_Modulus);
    fmpz_powm_ui(ghalf_, Int_Gen, step_size / 2, Int_Modulus);
    fmpz_set_ui(roots[0], 1);
    fmpz_set_ui(invroots[0], 1);
    fmpz_set_ui(roots2[0], 1);

    for (int i = 1; i < N; i++) {
        fmpz_mul(roots[i],roots[i-1],g_);
        fmpz_mul(invroots[i],invroots[i-1],ginv_);

        fmpz_mod(roots[i],roots[i],Int_Modulus);
        fmpz_mod(invroots[i],invroots[i],Int_Modulus);
    }

    for (int i = 1; i < 2 * N; i++) {
        fmpz_mul(roots2[i], roots2[i-1], ghalf_);
        fmpz_mod(roots2[i], roots2[i], Int_Modulus);
    }

    fmpz_clear(g_);
    fmpz_clear(ginv_);
    fmpz_clear(ghalf_);
}