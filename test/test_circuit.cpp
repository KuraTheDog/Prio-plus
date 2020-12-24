/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Demo FLINT program for incremental multimodular reduction and
    reconstruction using the Chinese Remainder Theorem.

    Example execution, done in test/ directory.

    gcc -c -std=c99 -o fft.o ../poly/fft.c && \
    g++ -c -std=c++11 -o test_circuit.o test_circuit.cpp && \
    g++ -c -std=c++11 -o fmpz_utils.o ../fmpz_utils.cpp && \
    g++ -c -std=c++11 -o share.o ../share.cpp && \
    g++ -o test_circuit fft.o test_circuit.o fmpz_utils.o share.o -lgmp -lflint -g && \
    ./test_circuit
*/

#include <iostream>
#include <gmpxx.h>

extern "C" {
  #include "flint/flint.h"
  #include "flint/fmpz.h"
  #include "flint/ulong_extras.h"
};

#include "../circuit.h"
#include "../share.h"
#include "../client.h"

int main(int argc, char* argv[])
{
  init_constants();
  Circuit* var_circuit = CheckVar();

  fmpz_t inp[2];
  fmpz_init(inp[0]);
  fmpz_init(inp[1]);

  fmpz_set_si(inp[0],9);
  fmpz_set_si(inp[1],81);

  std::cout << var_circuit->Eval(inp) << std::endl;
  ClientPacket p0, p1;
  
  share_polynomials(var_circuit,p0,p1);

  clear_constants();
  return 0;
}
