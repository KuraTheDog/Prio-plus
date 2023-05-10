#include "fmpz_utils.h"

#include <gmpxx.h>

extern "C" {
  #include "flint/flint.h"
  #include "flint/fmpz.h"
};

void init_set_fmpz_readonly(fmpz_t x, const fmpz_t in) {
  // // NOTE: This version fails on larger values.
  // mpz_t z; mpz_init(z);
  // fmpz_get_mpz(z, in);
  // fmpz_init_set_readonly(x, z);
  // mpz_clear(z);

  // Based on documentation.
  fmpz_t tmp; fmpz_init_set(tmp, in);  // Need to lose const
  fmpz_init_set_readonly(x, _fmpz_promote_val(tmp));
  _fmpz_demote_val(tmp);
  fmpz_clear(tmp);
}

void new_fmpz_array(fmpz_t** const arr, const size_t N) {
  fmpz_t* const out = (fmpz_t*) malloc(N * sizeof(fmpz_t));
  for (unsigned int i = 0; i < N; i++)
    fmpz_init_set_ui(out[i], 0);
  *arr = out;
}

void clear_fmpz_array(fmpz_t* arr, const size_t N) {
  for (unsigned int i = 0; i < N; i++)
    fmpz_clear(arr[i]);
  free(arr);
  arr = nullptr;
}

void copy_fmpz_array(fmpz_t* const dest, const fmpz_t* const src, const size_t N) {
  for (unsigned int i = 0; i < N; i++)
    fmpz_set(dest[i],src[i]);
}

// Turn bool (bit) array into fmpz_t
void fmpz_from_bool_array(fmpz_t x, const bool* const arr, const size_t n) {
  fmpz_zero(x);
  for (unsigned int i = 0; i < n; i++) {
    if (arr[i])
      fmpz_setbit(x, i);
  }
}

void fmpz_from_block(fmpz_t x, const emp::block &b, const size_t n) {
  fmpz_zero(x);
  for (unsigned int i = 0; i < n; i++)
    if ((1ULL << i) & (*((uint64_t*)&b)))
      fmpz_setbit(x, i);
}

void fmpz_mod_mul_si(fmpz_t a, const fmpz_t b, const slong c, const fmpz_mod_ctx_t ctx) {
  fmpz_mul_si(a, b, c);
  fmpz_mod(a, a, ctx->n);
}

void fmpz_mod_addmul(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx) {
  fmpz_t t; fmpz_init(t);
  fmpz_mod_mul(t, b, c, ctx);
  fmpz_mod_add(a, a, t, ctx);
  fmpz_clear(t);
}

void fmpz_mod_submul(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx) {
  fmpz_t t; fmpz_init(t);
  fmpz_mod_mul(t, b, c, ctx);
  fmpz_mod_sub(a, a, t, ctx);
  fmpz_clear(t);
}

void fmpz_mod_addmul_ui(fmpz_t a, const fmpz_t b, const ulong c, const fmpz_mod_ctx_t ctx) {
  fmpz_t t; fmpz_init(t);
  fmpz_mod_mul_ui(t, b, c, ctx);
  fmpz_mod_add(a, a, t, ctx);
  fmpz_clear(t);
}

void fmpz_mod_submul_ui(fmpz_t a, const fmpz_t b, const ulong c, const fmpz_mod_ctx_t ctx) {
  fmpz_t t; fmpz_init(t);
  fmpz_mod_mul_ui(t, b, c, ctx);
  fmpz_mod_sub(a, a, t, ctx);
  fmpz_clear(t);
}

void to_fsigned(fmpz_t x, const fmpz_t M) {
  fmpz_t half; fmpz_init(half);
  fmpz_cdiv_q_ui(half, M, 2);
  if (fmpz_cmp(x, half) > 0) {  // > N/2, negative
    fmpz_sub(x, x, M);
  }
  fmpz_clear(half);
}

int64_t get_fsigned(const fmpz_t x, const fmpz_t M) {
  fmpz_t tmp; fmpz_init_set(tmp, x);
  to_fsigned(tmp, M);
  int64_t ans = fmpz_get_si(tmp);
  fmpz_clear(tmp);
  return ans;
}
