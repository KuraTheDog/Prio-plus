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

#include <iostream>
#include <gmpxx.h>

extern "C" {
  #include "flint/flint.h"
  #include "flint/fmpz.h"
};

#include "../circuit.h"
#include "../share.h"
#include "../client.h"
#include "../server.h"

void test_CheckVar() {
  std::cout << "Testing CheckVar Eval and share_polynomials" << std::endl;
  Circuit* var_circuit = CheckVar();  // x^2 == y

  fmpz_t inp[2];
  fmpz_init(inp[0]);
  fmpz_init(inp[1]);

  fmpz_set_si(inp[0], 9);  // Set first input to 9
  fmpz_set_si(inp[1], 81);  // Set second to 81

  // Eval to verify 9^2 = 81
  /*
  Gate 0: input gate, x = 9
  Gate 1: input gate, y = 81
  Gate 2: mult gate, x^2 = 81
  Gate 3: negation gate, "-y" mod Int_Modulus
  Gate 4: addition gate, x^2 - y = 0, checked to be 0
  */
  bool eval = var_circuit->Eval(inp);
  std::cout << "Eval: " << eval << std::endl;

  /*
  One mult gate (#2). Next power of 2 is 2 (N).
  */
  ClientPacket p0, p1;
  share_polynomials(var_circuit,p0,p1);

  std::cout << "p0" << std::endl;
  p0->print();

  std::cout << "p1" << std::endl;
  p1->print();

  std::cout << "Running through validity checks" << std::endl;

  Circuit* var_circuit0 = CheckVar();
  Circuit* var_circuit1 = CheckVar();
  Checker* checker_0 = new Checker(var_circuit0, 0);
  Checker* checker_1 = new Checker(var_circuit1, 1);

  std::cout << "checker1 N: " << checker_1->N << std::endl;

  checker_0->setReq(p0);
  checker_1->setReq(p1);

  CheckerPreComp* pre0 = new CheckerPreComp(var_circuit);
  pre0->setCheckerPrecomp(inp[0]);  // TODO: what value? random?
  checker_0->evalPoly(pre0, p0);

  std::cout << "Checker 0:" << std::endl;
  std::cout << "  evalF = "; fmpz_print(checker_0->evalF); std::cout << std::endl;
  std::cout << "  evalG = "; fmpz_print(checker_0->evalG); std::cout << std::endl;
  std::cout << "  evalH = "; fmpz_print(checker_0->evalH); std::cout << std::endl;

  CheckerPreComp* pre1 = new CheckerPreComp(var_circuit);
  pre1->setCheckerPrecomp(inp[0]);  // TODO: what value? random?
  checker_1->evalPoly(pre1, p1);

  std::cout << "Checker 1:" << std::endl;
  std::cout << "  evalF = "; fmpz_print(checker_1->evalF); std::cout << std::endl;
  std::cout << "  evalG = "; fmpz_print(checker_1->evalG); std::cout << std::endl;
  std::cout << "  evalH = "; fmpz_print(checker_1->evalH); std::cout << std::endl;
}

int main(int argc, char* argv[])
{
  init_constants();

  test_CheckVar();

  clear_constants();
  return 0;
}
