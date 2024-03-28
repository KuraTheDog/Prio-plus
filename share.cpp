#include "share.h"

#include <flint/flint.h>
#include <flint/fmpz.h>

#include "constants.h"
#include "fmpz_utils.h"

unsigned int NextPowerOfTwo(const unsigned int n) {
    unsigned int ans = 1;
    while(n + 1 > ans)
        ans *= 2;
    return ans;
}

void SplitShare(const fmpz_t val, fmpz_t A, fmpz_t B) {
    fmpz_randm(A, seed, Int_Modulus);
    fmpz_sub(B, val, A);
    fmpz_mod(B, B, Int_Modulus);
}

Cor::Cor(const CorShare* const x, const CorShare* const y) : Cor() {
    fmpz_add(D, x->shareD, y->shareD);
    fmpz_mod(D, D, Int_Modulus);

    fmpz_add(E, x->shareE, y->shareE);
    fmpz_mod(E, E, Int_Modulus);
}

// Unused?
/*
void SplitShare(const fmpz_t val, fmpz_t A, fmpz_t B, const int num_bits) {
    // num_bits < 32
    uint64_t mod = 1L << num_bits;

    fmpz_t max_val;
    fmpz_init(max_val);
    fmpz_set_ui(max_val, mod);
    fmpz_randm(A, seed, max_val);
    fmpz_xor(B, A, val);
    fmpz_clear(max_val);
}
*/

BeaverTriple* NewBeaverTriple() {
    BeaverTriple* out = new BeaverTriple();

    fmpz_randm(out->A, seed, Int_Modulus);
    fmpz_randm(out->B, seed, Int_Modulus);
    fmpz_mul(out->C, out->A, out->B);
    fmpz_mod(out->C, out->C, Int_Modulus);

    return out;
}

void BeaverTripleShares(const BeaverTriple* const inp,
                        BeaverTripleShare* const out0,
                        BeaverTripleShare* const out1) {
    SplitShare(inp->A, out0->shareA, out1->shareA);
    SplitShare(inp->B, out0->shareB, out1->shareB);
    SplitShare(inp->C, out0->shareC, out1->shareC);
}

void makeLocalDaBit(DaBit* const bit0, DaBit* const bit1) {
    fmpz_t two;
    fmpz_init_set_si(two, 2);

    // random bit b
    fmpz_t bit;
    fmpz_init(bit);
    fmpz_randm(bit, seed, two);

    // random r
    /*
    If b >= r, then b - r doesn't roll over. This breaks assumptions for making b
    So we need r > b, so we randomly pick r between b+1 and modulus
    or, 0 and modulus - (b+1), then add b+1
    This case only happens odds ~1/p anyways, but since it's local we are free to adjust
    */
    fmpz_t adjust_mod; fmpz_init_set(adjust_mod, Int_Modulus);
    fmpz_t adjust; fmpz_init_set(adjust, bit); fmpz_add_si(adjust, adjust, 1);
    fmpz_sub(adjust_mod, adjust_mod, adjust);
    fmpz_randm(bit0->bp, seed, adjust_mod);
    fmpz_add(bit0->bp, bit0->bp, adjust);

    // b - r
    fmpz_sub(bit1->bp, bit, bit0->bp);
    fmpz_mod(bit1->bp, bit1->bp, Int_Modulus);

    // r mod 2. true if odd.
    bit0->b2 = fmpz_is_odd(bit0->bp);

    // 1 xor (b - r) mod 2 = (1 + b - r) mod 2, true if (b-r) even
    bit1->b2 = fmpz_is_even(bit1->bp);

    fmpz_clear(two);
    fmpz_clear(bit);
}

// Deprecated
void makeLocalEdaBit(EdaBit* const ebit0, EdaBit* const ebit1, const size_t n) {
    DaBit* bit0 = new DaBit();
    DaBit* bit1 = new DaBit();

    fmpz_zero(ebit0->r);
    fmpz_zero(ebit1->r);

    fmpz_t pow;
    fmpz_init_set_si(pow, 1);  // 2^i

    for (unsigned int i = 0; i < n; i++) {
        makeLocalDaBit(bit0, bit1);

        ebit0->b[i] = bit0->b2;
        ebit1->b[i] = bit1->b2;

        fmpz_addmul(ebit0->r, pow, bit0->bp);
        fmpz_mod(ebit0->r, ebit0->r, Int_Modulus);
        fmpz_addmul(ebit1->r, pow, bit1->bp);
        fmpz_mod(ebit1->r, ebit1->r, Int_Modulus);

        fmpz_mul_si(pow, pow, 2);
        fmpz_mod(pow, pow, Int_Modulus);
    }
    fmpz_clear(pow);
    delete bit0;
    delete bit1;
}
