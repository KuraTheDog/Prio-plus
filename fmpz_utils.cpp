#include "fmpz_utils.h"

#include <gmpxx.h>

extern "C" {
  #include "flint/flint.h"
  #include "flint/fmpz.h"
};

void init_fmpz_array(fmpz_t *arr, int N){
    for(int i = 0; i < N; i++)
        fmpz_init(arr[i]);
}

void set_zero_fmpz_array(fmpz_t* arr, int N){
    for(int i = 0; i < N; i++)
        fmpz_set_ui(arr[i],0);
}

void new_fmpz_array(fmpz_t** arr, int N){
    fmpz_t* out = (fmpz_t*) malloc(N * sizeof(fmpz_t));
    init_fmpz_array(out,N);
    set_zero_fmpz_array(out,N);
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
