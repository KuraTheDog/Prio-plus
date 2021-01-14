#include "../share.h"
#include <assert.h>

void testSplitShare() {
    fmpz_t num, sum;
    init_constants();

    fmpz_init_set_ui(num,1000);
    fmpz_init(sum);

    fmpz_t* ans = SplitShare(num);

    fmpz_add(sum,ans[0],ans[1]);
    fmpz_mod(sum,sum,Int_Modulus);

    assert(fmpz_equal(num,sum) == 1);

    fmpz_clear(num);
    fmpz_clear(sum);
    fmpz_clear(ans[0]);
    fmpz_clear(ans[1]);
    free(ans);
}

void testBeaverTriple() {
    BeaverTriple* t = NewBeaverTriple();
    BeaverTripleShare* out = BeaverTripleShares(t);
    fmpz_t A, B, C, prod;
    fmpz_init(A);
    fmpz_init(B);
    fmpz_init(C);
    fmpz_init(prod);

    fmpz_add(A,out[0].shareA,out[1].shareA);
    fmpz_mod(A,A,Int_Modulus);
    fmpz_add(B,out[0].shareB,out[1].shareB);
    fmpz_mod(B,B,Int_Modulus);
    fmpz_add(C,out[0].shareC,out[1].shareC);
    fmpz_mod(C,C,Int_Modulus);
    fmpz_mul(prod,A,B);
    fmpz_mod(prod,prod,Int_Modulus);

    assert(fmpz_equal(A,t->A) == 1 and fmpz_equal(B,t->B) == 1 and fmpz_equal(C,t->C) == 1 and fmpz_equal(prod,t->C) == 1);

    fmpz_clear(A);
    fmpz_clear(B);
    fmpz_clear(C);
    delete t;
    delete[] out;
}

int main() {
    testSplitShare();
    testBeaverTriple();


    return 0;
}