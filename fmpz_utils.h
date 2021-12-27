#ifndef FMPZ_UTILS_H
#define FMPZ_UTILS_H

#include <emp-tool/emp-tool.h>
#include <gmpxx.h>

extern "C" {
  #include "flint/flint.h"
  #include "flint/fmpz.h"
};

// Note: also array is initialized to 0.
void new_fmpz_array(fmpz_t** arr, const size_t N);

void clear_fmpz_array(fmpz_t* arr, const size_t N);

void copy_fmpz_array(fmpz_t* dest, const fmpz_t* const src, const size_t N);

void fmpz_from_bool_array(fmpz_t x, const bool* const arr, const size_t n);

// n bit fmpz, with n up to 256
void fmpz_from_block(fmpz_t x, const emp::block &b, const size_t n);

#endif