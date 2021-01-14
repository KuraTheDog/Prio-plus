#ifndef _FFT__H
#define _FFT__H

#include <stdbool.h>

#include <flint/fmpz_mod_poly.h>

fmpz_t *fft_interpolate(const fmpz_t mod, const int nPoints,
    const fmpz_t *roots, const fmpz_t *ys, const bool invert);

#endif
