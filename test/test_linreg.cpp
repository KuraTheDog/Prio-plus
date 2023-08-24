#undef NDEBUG
#include <assert.h>
#include <iostream>

#include "../circuit.h"
#include "../client.h"
#include "../server.h"
#include "../share.h"

void test_CheckLinReg() {
  std::cout << "Testing CheckLinReg Eval and share_polynomials" << std::endl;
  Circuit* linreg_circuit = CheckLinReg(2);

  // Set inputs
  fmpz_t inp[4];
  fmpz_init(inp[0]); fmpz_set_si(inp[0], 2);  // x
  fmpz_init(inp[1]); fmpz_set_si(inp[1], 3);  // y
  fmpz_init(inp[2]); fmpz_set_si(inp[2], 4);  // x^2
  fmpz_init(inp[3]); fmpz_set_si(inp[3], 6);  // xy

  // Check x * y = xy, x * x = x^2
  bool eval = linreg_circuit->Eval(inp);
  std::cout << "Eval: " << eval << std::endl;

  // Setup packets
  size_t N = NextPowerOfTwo(linreg_circuit->NumMulGates());
  ClientPacket* p0 = new ClientPacket(linreg_circuit->NumMulGates());
  ClientPacket* p1 = new ClientPacket(linreg_circuit->NumMulGates());
  share_polynomials(linreg_circuit, p0, p1);
  delete linreg_circuit;
  std::cout << "p0" << std::endl; p0->print();
  std::cout << "p1" << std::endl; p1->print();

  // Set arithmetic shares
  fmpz_t inp0[4];
  fmpz_t inp1[4];
  fmpz_init(inp0[0]); fmpz_init(inp1[0]);
  SplitShare(inp[0], inp0[0], inp1[0]);
  fmpz_init(inp0[1]); fmpz_init(inp1[1]);
  SplitShare(inp[1], inp0[1], inp1[1]);
  fmpz_init(inp0[2]); fmpz_init(inp1[2]);
  SplitShare(inp[2], inp0[2], inp1[2]);
  fmpz_init(inp0[3]); fmpz_init(inp1[3]);
  SplitShare(inp[3], inp0[3], inp1[3]);

  std::cout << "------ Running through validity checks" << std::endl;

  fmpz_t randomX;
  fmpz_init(randomX);
  fmpz_randm(randomX, seed, Int_Modulus);
  std::cout << "Random X: "; fmpz_print(randomX); std::cout << std::endl;

  Circuit* linreg_circuit0 = CheckLinReg(2);
  MultCheckPreComp* pre0 = new MultCheckPreComp(N);
  pre0->setEvalPoint(randomX);

  Checker* checker_0 = new Checker(linreg_circuit0, 0, p0, pre0, inp0, true);

  std::cout << "-=-=-=-=-=-" << std::endl;

  Circuit* linreg_circuit1 = CheckLinReg(2);
  MultCheckPreComp* pre1 = new MultCheckPreComp(N);
  pre1->setEvalPoint(randomX);

  Checker* checker_1 = new Checker(linreg_circuit1, 1, p1, pre1, inp1, true);

  auto corshare0 = checker_0->CorFn();
  auto corshare1 = checker_1->CorFn();

  Cor* cor0 = new Cor(corshare0, corshare1);
  Cor* cor1 = new Cor(corshare0, corshare1);

  fmpz_t out0, out1;
  fmpz_init(out0);
  fmpz_init(out1);

  checker_0->OutShare(out0, cor0);
  checker_1->OutShare(out1, cor1);

  bool result = AddToZero(out0, out1);

  std::cout << "out0 : "; fmpz_print(out0); std::cout << ", out1 : "; fmpz_print(out1); std::cout << std::endl;

  assert(result == 1);
  std::cout << "Result : " << std::boolalpha << result << std::endl;

  std::cout << "^v^v^ Shared validation: " << std::endl;
  fmpz_t tmp, rgr;
  fmpz_init(rgr);
  fmpz_init(tmp);
  fmpz_mod_add(tmp, checker_0->evalF, checker_1->evalF, mod_ctx);
  std::cout << "f(r) = "; fmpz_print(tmp); std::cout << std::endl;
  fmpz_mod_add(rgr, checker_0->evalG, checker_1->evalG, mod_ctx);
  std::cout << "r * g(r) = "; fmpz_print(rgr); std::cout << std::endl;
  fmpz_mod_mul(tmp, tmp, rgr, mod_ctx);
  std::cout << "r * f(r) * g(r) = "; fmpz_print(tmp); std::cout << std::endl;

  fmpz_mod_add(tmp, checker_0->evalH, checker_1->evalH, mod_ctx);
  std::cout << "r * h(r) = "; fmpz_print(tmp); std::cout << std::endl;

  fmpz_clear(out0);
  fmpz_clear(out1);
  delete cor0;
  delete cor1;
  delete corshare0;
  delete corshare1;
  delete pre1;
  delete checker_1;
  delete linreg_circuit1;
  delete pre0;
  delete checker_0;
  delete linreg_circuit0;
  fmpz_clear(randomX);

  delete p0;
  delete p1;
  for (unsigned int i = 0; i < 4; i++) {
    fmpz_clear(inp[i]);
    fmpz_clear(inp0[i]);
    fmpz_clear(inp1[i]);
  }
}

int main(int argc, char* argv[]) {
  init_constants();

  test_CheckLinReg();

  RootManager(1).clearCache();
  clear_constants();
  return 0;
}
