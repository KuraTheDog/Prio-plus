#ifndef FMPZ_UTILS_H
#define FMPZ_UTILS_H

#include <emp-tool/emp-tool.h>
#include <gmpxx.h>

extern "C" {
  #include "flint/flint.h"
  #include "flint/fmpz.h"
  #include "flint/fmpz_mod.h"
};

void init_set_fmpz_readonly(fmpz_t x, const fmpz_t in);

void new_fmpz_array(fmpz_t** arr, const size_t N);

void clear_fmpz_array(fmpz_t* arr, const size_t N);

void copy_fmpz_array(fmpz_t* dest, const fmpz_t* const src, const size_t N);

void fmpz_from_bool_array(fmpz_t x, const bool* const arr, const size_t n);

// n bit fmpz, with n up to 256
void fmpz_from_block(fmpz_t x, const emp::block &b, const size_t n);

// Not currently part of "fmpz_mod" library.
// Convenience functions
void fmpz_mod_mul_si(fmpz_t a, const fmpz_t b, const slong c, const fmpz_mod_ctx_t ctx);
void fmpz_mod_addmul(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx);
void fmpz_mod_addmul_ui(fmpz_t a, const fmpz_t b, const ulong c, const fmpz_mod_ctx_t ctx);

#endif