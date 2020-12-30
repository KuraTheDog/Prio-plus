#include "fmpz_utils.h"

#include <gmpxx.h>

extern "C" {
  #include "flint/flint.h"
  #include "flint/fmpz.h"
};

void new_fmpz_array(fmpz_t** arr, int N){
    fmpz_t* out = (fmpz_t*) malloc(N * sizeof(fmpz_t));
    for (int i = 0; i < N; i++)
        fmpz_init_set_ui(out[i], 0);
    *arr = out;
}

void clear_fmpz_array(fmpz_t* arr, int N){
    for(int i = 0; i < N; i++)
        fmpz_clear(arr[i]);
    free(arr);
}

void copy_fmpz_array(fmpz_t* dest, fmpz_t* src, int N){
    for(int i = 0; i < N; i++)
        fmpz_set(dest[i],src[i]);
}
