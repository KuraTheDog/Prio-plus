#include "fmpz_utils.h"

#include <gmpxx.h>

extern "C" {
  #include "flint/flint.h"
  #include "flint/fmpz.h"
};

void new_fmpz_array(fmpz_t** arr, const int N) {
    fmpz_t* out = (fmpz_t*) malloc(N * sizeof(fmpz_t));
    for (int i = 0; i < N; i++)
        fmpz_init_set_ui(out[i], 0);
    *arr = out;
}

void clear_fmpz_array(fmpz_t* arr, const int N) {
    for (int i = 0; i < N; i++)
        fmpz_clear(arr[i]);
    free(arr);
}

void copy_fmpz_array(fmpz_t* dest, const fmpz_t* src, const int N) {
    for (int i = 0; i < N; i++)
        fmpz_set(dest[i],src[i]);
}

bool get_fmpz_bit(const fmpz_t x, const size_t n) {
    fmpz_t mask;
    fmpz_init_set_ui(mask, 1 << n);
    fmpz_and(mask, x, mask);  // just nth bit of x
    bool ans = (1 - fmpz_is_zero(mask));  // 0 if mask is zero
    fmpz_clear(mask);
    return ans;
}
