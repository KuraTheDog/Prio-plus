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

    Except not really
*/
/*
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
  Circuit* var_circuit = CheckVar();  // x^2 == y

  fmpz_t inp[2];
  fmpz_init(inp[0]);
  fmpz_init(inp[1]);

  fmpz_set_si(inp[0], 9);  // Set first input to 9
  fmpz_set_si(inp[1], 81);  // Set second to 81

  // Eval to verify 9^2 = 81
  /*
  Gate 0: input 1 = x = 9
  Gate 1: input 2 = y = 81
  Gate 2: x^2 = 81
  Gate 3: "-y" mod Int_Modulus
  Gate 4: x^2 - y = 0, checked to be 0
  */
  bool eval = var_circuit->Eval(inp);
  std::cout << "Eval: " << eval << std::endl;

  ClientPacket p0, p1;
  share_polynomials(var_circuit,p0,p1);

  std::cout << "p0" << std::endl;
  p0->print();

  std::cout << "p1" << std::endl;
  p1->print();

  clear_constants();
  return 0;
}
