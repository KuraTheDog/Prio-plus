#include "share.h"

#include "constants.h"
#include "fmpz_utils.h"

unsigned int NextPowerOfTwo(const unsigned int n) {
    unsigned int ans = 1;
    while(n + 1 > ans)
        ans *= 2;
    return ans;
}

size_t accumulate(const size_t num_inputs, const size_t num_values,
                  const fmpz_t* const values, const bool* const valid,
                  fmpz_t* const ans) {
  size_t num_valid = 0;

  for (unsigned int i = 0; i < num_inputs; i++) {
    if (!valid[i])
      continue;
    for (unsigned int j = 0; j < num_values; j++) {
      fmpz_mod_add(ans[j], ans[j], values[i * num_values + j], mod_ctx);
    }
    num_valid++;
  }
  return num_valid;
}

void SplitShare(const fmpz_t val, fmpz_t A, fmpz_t B) {
    fmpz_randm(A, seed, Int_Modulus);
    fmpz_mod_sub(B, val, A, mod_ctx);
}

Cor::Cor(const Cor* const x, const Cor* const y) : Cor() {
    fmpz_mod_add(D, x->D, y->D, mod_ctx);
    fmpz_mod_add(E, x->E, y->E, mod_ctx);
}

void makeLocalBoolTriple(BooleanBeaverTriple* const out0, BooleanBeaverTriple* const out1) {
    fmpz_t r; fmpz_init(r);
    fmpz_randbits(r, seed, 5);
    out0->a = fmpz_tstbit(r, 0);
    out1->a = fmpz_tstbit(r, 1);
    out0->b = fmpz_tstbit(r, 2);
    out1->b = fmpz_tstbit(r, 3);
    out0->c = fmpz_tstbit(r, 4);
    fmpz_clear(r);

    out1->c = out0->c ^ ((out0->a ^ out1->a) & (out0->b ^ out1->b));
}

void makeLocalTriple(BeaverTriple* const out0, BeaverTriple* const out1) {
    fmpz_randm(out0->A, seed, Int_Modulus);
    fmpz_randm(out1->A, seed, Int_Modulus);
    fmpz_randm(out0->B, seed, Int_Modulus);
    fmpz_randm(out1->B, seed, Int_Modulus);
    fmpz_randm(out0->C, seed, Int_Modulus);
    // (a + a)(b + b) = c + c
    fmpz_t tmp; fmpz_init(tmp);
    fmpz_mod_add(out1->C, out0->A, out1->A, mod_ctx);
    fmpz_mod_add(tmp, out0->B, out1->B, mod_ctx);
    fmpz_mod_mul(out1->C, out1->C, tmp, mod_ctx);
    fmpz_mod_sub(out1->C, out1->C, out0->C, mod_ctx);
}

void makeLocalAltTriple(AltTriple* const out0, AltTriple* const out1) {
    fmpz_randm(out0->AB, seed, Int_Modulus);
    fmpz_randm(out1->AB, seed, Int_Modulus);
    fmpz_randm(out0->C, seed, Int_Modulus);
    // a * b = c
    fmpz_mod_mul(out1->C, out0->AB, out1->AB, mod_ctx);
    fmpz_mod_sub(out1->C, out1->C, out0->C, mod_ctx);
}

void makeLocalDaBit(DaBit* const bit0, DaBit* const bit1) {
    // random bit b
    fmpz_t bit; fmpz_init(bit);
    fmpz_randbits(bit, seed, 1);

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
    fmpz_mod_add(bit0->bp, bit0->bp, adjust, mod_ctx);

    // b - r
    fmpz_mod_sub(bit1->bp, bit, bit0->bp, mod_ctx);

    // r mod 2. true if odd.
    bit0->b2 = fmpz_is_odd(bit0->bp);

    // 1 xor (b - r) mod 2 = (1 + b - r) mod 2, true if (b-r) even
    bit1->b2 = fmpz_is_even(bit1->bp);

    fmpz_clear(bit);
}

// Deprecated
void makeLocalEdaBit(EdaBit* const ebit0, EdaBit* const ebit1, const size_t n) {
    DaBit* const bit0 = new DaBit();
    DaBit* const bit1 = new DaBit();

    fmpz_zero(ebit0->r);
    fmpz_zero(ebit1->r);

    fmpz_t pow;
    fmpz_init_set_si(pow, 1);  // 2^i

    for (unsigned int i = 0; i < n; i++) {
        makeLocalDaBit(bit0, bit1);

        ebit0->b[i] = bit0->b2;
        ebit1->b[i] = bit1->b2;

        fmpz_mod_addmul(ebit0->r, pow, bit0->bp, mod_ctx);
        fmpz_mod_addmul(ebit1->r, pow, bit1->bp, mod_ctx);

        fmpz_mod_mul_si(pow, pow, 2, mod_ctx);
    }
    fmpz_clear(pow);
    delete bit0;
    delete bit1;
}
