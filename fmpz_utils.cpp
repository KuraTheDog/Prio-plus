#include "fmpz_utils.h"

#include <gmpxx.h>

extern "C" {
  #include "flint/flint.h"
  #include "flint/fmpz.h"
};

void init_set_fmpz_readonly(fmpz_t x, const fmpz_t in) {
  mpz_t z; mpz_init(z);
  fmpz_get_mpz(z, in);
  fmpz_init_set_readonly(x, z);
  mpz_clear(z);

  // Maybe more efficient, loses const. We don't do this enough to care.
  // Based on documentation.
  // fmpz_init_set_readonly(x, _fmpz_promote_val(in));
  // _fmpz_demote_val(in);
}

void new_fmpz_array(fmpz_t** const arr, const size_t N) {
  fmpz_t* const out = (fmpz_t*) malloc(N * sizeof(fmpz_t));
  for (unsigned int i = 0; i < N; i++)
    fmpz_init_set_ui(out[i], 0);
  *arr = out;
}

void clear_fmpz_array(fmpz_t* const arr, const size_t N) {
  for (unsigned int i = 0; i < N; i++)
    fmpz_clear(arr[i]);
  free(arr);
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
