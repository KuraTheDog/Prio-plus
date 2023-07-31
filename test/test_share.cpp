#undef NDEBUG
#include <assert.h>

#include "../constants.h"
#include "../share.h"

void testSplitShare() {
    fmpz_t num; fmpz_init_set_ui(num, 1000);
    fmpz_t a; fmpz_init(a);
    fmpz_t b; fmpz_init(b);

    SplitShare(num, a, b);

    fmpz_t sum; fmpz_init(sum);
    fmpz_mod_add(sum, a, b, mod_ctx);

    assert(fmpz_equal(num, sum) == 1);

    fmpz_clear(num);
    fmpz_clear(sum);
    fmpz_clear(a);
    fmpz_clear(b);
}

void testBeaverTriple() {
    BeaverTriple* t0 = new BeaverTriple();
    BeaverTriple* t1 = new BeaverTriple();
    makeLocalTriple(t0, t1);

    fmpz_t A, B, C, prod;
    fmpz_init(A);
    fmpz_init(B);
    fmpz_init(C);
    fmpz_init(prod);

    fmpz_mod_add(A, t0->A, t1->A, mod_ctx);
    fmpz_mod_add(B, t0->B, t1->B, mod_ctx);
    fmpz_mod_add(C, t0->C, t1->C, mod_ctx);
    fmpz_mod_mul(prod, A, B, mod_ctx);

    assert(fmpz_equal(prod, C) == 1);

    fmpz_clear(A);
    fmpz_clear(B);
    fmpz_clear(C);
    delete t0;
    delete t1;
}

int main() {
    init_constants();

    testSplitShare();
    testBeaverTriple();

    return 0;
}
