#ifndef _FFT__H
#define _FFT__H

#include <stdbool.h>

#include <flint/fmpz.h>

fmpz_t* fft_interpolate(const fmpz_t mod, const int nPoints,
    const fmpz_t* const roots, const fmpz_t* const ys, const bool invert);

#endif
