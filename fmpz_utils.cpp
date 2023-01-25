#include "fmpz_utils.h"

#include <gmpxx.h>

extern "C" {
  #include "flint/flint.h"
  #include "flint/fmpz.h"
};

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
