#undef NDEBUG
#include <assert.h>
#include <iostream>

#include "../constants.h"
#include "../interp.h"
#include "../fmpz_utils.h"

void test_roots(const size_t N) {
  std::cout << "Testing roots: " << N << std::endl;
  RootManager manager = RootManager(N);

  fmpz_t* roots = manager.getRoots();
  fmpz_t* roots2 = manager.getRoots2();
  fmpz_t* rootsInv = manager.getRootsInv();

  // all start with 1
  assert(fmpz_equal_ui(roots[0], 1));
  assert(fmpz_equal_ui(roots2[0], 1));
  assert(fmpz_equal_ui(rootsInv[0], 1));

  assert(fmpz_equal(roots2[2], roots[1]));
  fmpz_t tmp; fmpz_init(tmp);
  fmpz_mul(tmp, roots[1], rootsInv[1]);
  fmpz_mod(tmp, tmp, Int_Modulus);
  assert(fmpz_equal_ui(tmp, 1));
  // deal with unused
  fmpz_set(tmp, roots2[0]);

  fmpz_clear(tmp);
  manager.clearCache();
}

void test_basic_fft(const size_t N) {
  std::cout << "Testing basic fft: " << N << std::endl;
  RootManager manager = RootManager(N);

  fmpz_t* points; new_fmpz_array(&points, N);
  for (unsigned int i = 0; i < N; i++) {
    fmpz_randm(points[i], seed, Int_Modulus);
  }

  fmpz_t* coeffs = interpolate_inv(N, points);
  fmpz_t* evals = interpolate_N(N, coeffs);

  // poly though points, then eval at points. ensure matches.
  for (unsigned int i = 0; i < N; i++) {
    assert(fmpz_equal(points[i], evals[i]));
  }

  clear_fmpz_array(points, N);
  clear_fmpz_array(coeffs, N);
  clear_fmpz_array(evals, N);
}

void randomize(fmpz_t* const arr, const size_t N) {
  for (unsigned int i = 0; i < N; i++) {
    fmpz_randm(arr[i], seed, Int_Modulus);
  }
}

void eval(const fmpz_t* const coeffs, const fmpz_t x, const size_t N, fmpz_t out) {
  fmpz_t xpow; fmpz_init_set_ui(xpow, 1);
  fmpz_set_ui(out, 0);
  for (unsigned int i = 0; i < N; i++) {
    fmpz_addmul(out, coeffs[i], xpow);
    fmpz_mul(xpow, xpow, x);
    fmpz_mod(xpow, xpow, Int_Modulus);
    fmpz_mod(out, out, Int_Modulus);
  }
}

void test_precompute(const size_t N) {
  std::cout << "Testing MultCheckPreComp: " << N << std::endl;
  // Checks eval at 1 is first point.
  // More complex to check other points.

  MultCheckPreComp* chk = new MultCheckPreComp(N);
  fmpz_t tmp; fmpz_init(tmp);

  fmpz_t* points; new_fmpz_array(&points, N);

  fmpz_set_ui(tmp, 1);
  chk->setEvalPoint(tmp);

  // check evals at 1 matches first point
  for (unsigned int i = 0; i < 5; i++) {
    randomize(points, N);
    chk->Eval(points, tmp);
    assert(fmpz_equal(points[0], tmp));
  }

  // Actually check a poly
  // random coeffs, manually eval at points and at an eval point
  fmpz_t* coeffs; new_fmpz_array(&coeffs, N);
  randomize(coeffs, N);
  fmpz_set_ui(tmp, 42);
  chk->setEvalPoint(tmp);
  fmpz_t expected; fmpz_init(expected);
  const fmpz_t* const roots = RootManager(N).getRoots();
  eval(coeffs, tmp, N, expected);
  for (unsigned int i = 0; i < N; i++) {
    eval(coeffs, roots[i], N, points[i]);
  }
  chk->Eval(points, tmp);
  assert(fmpz_equal(tmp, expected));

  fmpz_clear(expected);
  fmpz_clear(tmp);
  clear_fmpz_array(points, N);
  clear_fmpz_array(coeffs, N);

  delete chk;
}

int main(int argc, char** argv){
  init_constants();

  test_roots(4);
  test_roots(8);

  test_basic_fft(4);
  test_basic_fft(8);

  test_precompute(4);
  test_precompute(8);

  RootManager(4).clearCache();

  clear_constants();
}