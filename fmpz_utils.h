#ifndef FMPZ_UTILS_H
#define FMPZ_UTILS_H

#include <gmpxx.h>

extern "C" {
  #include "flint/flint.h"
  #include "flint/fmpz.h"
};

void new_fmpz_array(fmpz_t** arr, const int N);

void clear_fmpz_array(fmpz_t* arr, const int N);

void copy_fmpz_array(fmpz_t* dest, const fmpz_t* src, const int N);

bool get_fmpz_bit(const fmpz_t x, const size_t n);

void fmpz_from_bool_array(fmpz_t x, const bool* arr, const size_t n);

#endif