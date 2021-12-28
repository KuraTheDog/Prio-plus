#include "fmpz_utils.h"

#include <gmpxx.h>

extern "C" {
  #include "flint/flint.h"
  #include "flint/fmpz.h"
};

void new_fmpz_array(fmpz_t** arr, const size_t N) {
    fmpz_t* out = (fmpz_t*) malloc(N * sizeof(fmpz_t));
    for (unsigned int i = 0; i < N; i++)
        fmpz_init_set_ui(out[i], 0);
    *arr = out;
}

void clear_fmpz_array(fmpz_t* arr, const size_t N) {
    for (unsigned int i = 0; i < N; i++)
        fmpz_clear(arr[i]);
    free(arr);
}

void copy_fmpz_array(fmpz_t* dest, const fmpz_t* const src, const size_t N) {
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

int64_t get_fsigned(const fmpz_t x, const fmpz_t M) {
  fmpz_t tmp; fmpz_init(tmp);
  fmpz_cdiv_q_ui(tmp, M, 2);
  int64_t ans;

  if (fmpz_cmp(x, tmp) > 0) {  // > N/2, negative
    fmpz_sub(tmp, x, M);
    ans = fmpz_get_si(tmp);
  } else {
    ans = fmpz_get_si(x);
  }
  fmpz_clear(tmp);
  return ans;
}
