#include "constants.h"

#include <iostream>

#include "fmpz_utils.h"

fmpz_t Int_Modulus;
fmpz_mod_ctx_t mod_ctx;
fmpz_t Int_Gen;
size_t mod_size;
flint_rand_t seed;
fmpz_t *roots = nullptr, *invroots = nullptr, *roots2 = nullptr;
size_t num_roots;

void init_constants() {
    fmpz_t tmp; fmpz_init(tmp);

    fmpz_set_str(tmp, Int_Modulus_str.c_str(), 16);
    init_set_fmpz_readonly(Int_Modulus, tmp);

    fmpz_mod_ctx_init(mod_ctx, Int_Modulus);

    mod_size = fmpz_size(Int_Modulus);

    fmpz_set_str(tmp, Int_Gen_str.c_str(), 16);
    init_set_fmpz_readonly(Int_Gen, tmp);

    fmpz_clear(tmp);

    std::cout << "Init constants:\n";
    std::cout << "  Int_Modulus = "; fmpz_print(Int_Modulus); std::cout << "\n";
    std::cout << "  Int_Gen = "; fmpz_print(Int_Gen); std::cout << "\n";
    std::cout << "  twoOrder = " << twoOrder << std::endl;

    flint_randinit(seed);
}

void clear_constants() {
    flint_randclear(seed);
    fmpz_clear_readonly(Int_Modulus);
    fmpz_mod_ctx_clear(mod_ctx);
    fmpz_clear_readonly(Int_Gen);
}

void init_roots(const size_t N) {
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

    const int step_size = (1 << twoOrder)/N;  // 2^(twoOrder - log_2 N)
    fmpz_t g_, ginv_, ghalf_;
    fmpz_init(g_);      // Generator (Int_Gen)
    fmpz_init(ginv_);   // Inverse of g_
    fmpz_init(ghalf_);  // g_^(step/2), for 2N roots.

    fmpz_mod_inv(ginv_, Int_Gen, mod_ctx);

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

    for (unsigned int i = 1; i < N; i++) {
        fmpz_mod_mul(roots[i], roots[i-1], g_, mod_ctx);
        fmpz_mod_mul(invroots[i], invroots[i-1], ginv_, mod_ctx);
    }

    for (unsigned int i = 1; i < 2 * N; i++) {
        fmpz_mod_mul(roots2[i], roots2[i-1], ghalf_, mod_ctx);
    }

    fmpz_clear(g_);
    fmpz_clear(ginv_);
    fmpz_clear(ghalf_);
}
