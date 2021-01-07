#ifndef _FFT__H
#define _FFT__H

#include <stdbool.h>

#include <flint/fmpz_mod_poly.h>

fmpz_t *fft_interpolate(fmpz_t mod, int nPoints, 
    fmpz_t *roots, fmpz_t *ys, bool invert);

#endif
