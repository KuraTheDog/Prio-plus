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
    fmpz_add(sum, a, b);
    fmpz_mod(sum, sum, Int_Modulus);

    assert(fmpz_equal(num, sum) == 1);

    fmpz_clear(num);
    fmpz_clear(sum);
    fmpz_clear(a);
    fmpz_clear(b);
}

void testBeaverTriple() {
    BeaverTriple* t0 = new BeaverTriple();
    BeaverTriple* t1 = new BeaverTriple();
    NewBeaverTriples(t0, t1);

    fmpz_t A, B, C, prod;
    fmpz_init(A);
    fmpz_init(B);
    fmpz_init(C);
    fmpz_init(prod);

    fmpz_add(A, t0->A, t1->A);
    fmpz_mod(A, A, Int_Modulus);
    fmpz_add(B, t0->B, t1->B);
    fmpz_mod(B, B, Int_Modulus);
    fmpz_add(C, t0->C, t1->C);
    fmpz_mod(C, C, Int_Modulus);
    fmpz_mul(prod, A, B);
    fmpz_mod(prod, prod, Int_Modulus);

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