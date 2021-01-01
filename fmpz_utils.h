#ifndef FMPZ_UTILS_H
#define FMPZ_UTILS_H

#include <gmpxx.h>

extern "C" {
  #include "flint/flint.h"
  #include "flint/fmpz.h"
};

void new_fmpz_array(fmpz_t** arr, int N);

void clear_fmpz_array(fmpz_t* arr, int N);

void copy_fmpz_array(fmpz_t* dest, fmpz_t* src, int N);

#endif