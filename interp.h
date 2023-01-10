#ifndef INTERP_H
#define INTERP_H

#include "unordered_map"

#include "constants.h"
#include "fmpz_utils.h"
#include "utils.h"

extern "C" {
    #include "poly/poly_batch.h"
    #include "poly/poly_once.h"
}

/*

Polynomial interpolation using FFT tricks, for identity testing.

Evaluation x points are fixed, and based on roots of unity.
Modulus is chosen for this, so that there's a 2^k generator.
So built to support various N = 2^k interpolation points.
RootManager caches them for different N

Then, the final evaluation point can be reused often.
Hence more precomputation can be done with it fixed, and occasionally updated
Eval then takes in the Y values corresponding to the X points, and it finishes
finding the corresponding poly, and returns it evaluated at the point.

MultCheckPreComp:
The common use case is mass testing multiplication.
Namely, for lists of {a}, {b}, {c}, ensure all a_i * b_i = c_i
For this, there is f(x_i) = a_i, g for b, and h for c.
Then identity test f * g = h at random x gives that all mults hold with high probability
For degrees, f and g are deg N, then h has to be deg 2N
Hence this gives an easy framework for having both N and 2N poly interp.
*/

/* Cache of roots of unity
  Nth roots of unity, used for FFT, where N is a power of 2
  Then we have roots = 1, r, r^2, ..., r^{N-1}, and their inverses.
  root2 are the 2Nth roots of unity.

  Use: RootManager(N).getRoots();
  computes roots around N if not cached, then fetches.
*/
class RootManager {
  static std::unordered_map<size_t, fmpz_t*> roots;
  static std::unordered_map<size_t, fmpz_t*> invroots;
  static std::unordered_map<size_t, fmpz_t*> roots2;
  static std::unordered_map<size_t, fmpz_t*> invroots2;

  const size_t N;

  void addRoots(const size_t N);

public:
  // Only powers of 2
  RootManager(const size_t N) : N(N) {
    if (not isPowerOfTwo(N))
      error_exit("RootManager only takes powers of 2");
    addRoots(N);
  };

  fmpz_t* getRoots() const {
    return roots.at(N);
  };

  fmpz_t* getRoots2() const {
    return roots2.at(N);
  };

  fmpz_t* getRootsInv() const {
    return invroots.at(N);
  };

  fmpz_t* getRoots2Inv() const {
    return invroots2.at(N);
  };

  void clearCache() {
    for (auto it = roots.begin(); it != roots.end(); ++it) {
      clear_fmpz_array(it->second, it->first);
    }
    roots.clear();
    for (auto it = invroots.begin(); it != invroots.end(); ++it) {
      clear_fmpz_array(it->second, it->first);
    }
    invroots.clear();
    for (auto it = roots2.begin(); it != roots2.end(); ++it) {
      clear_fmpz_array(it->second, it->first);
    }
    roots2.clear();
    for (auto it = invroots2.begin(); it != invroots2.end(); ++it) {
      clear_fmpz_array(it->second, it->first);
    }
    invroots2.clear();
  }

  ~RootManager() {};
};

/* "Normal" interpolate is just mass evaluate, with points=coeffs
Inverse interpolate finds coefficients.
*/
fmpz_t* interpolate_N(const size_t N, const fmpz_t* const points);
fmpz_t* interpolate_inv(const size_t N, const fmpz_t* const points);
fmpz_t* interpolate_2N(const size_t N, const fmpz_t* const points);
fmpz_t* interpolate_2N_inv(const size_t N, const fmpz_t* const points);


// Precompute on interpolation X values
struct BatchPre {
  precomp_t pre;
  BatchPre(const fmpz_t* const xPointsIn, const size_t n) {
    poly_batch_precomp_init(&pre, Int_Modulus, n, xPointsIn);
  }
  ~BatchPre() {
    poly_batch_precomp_clear(&pre);
  }
};

// Build precompute for evaluation point
struct PreX {
  const BatchPre* const batchPre;
  precomp_x_t pre;

  PreX(const BatchPre* const b, const fmpz_t x)
  : batchPre(b) {
    precomp_x_init(&pre, &batchPre->pre, x);
  }

  ~PreX() {
    precomp_x_clear(&pre);
  }

  void Eval(const fmpz_t* const yValues, fmpz_t out) {
    precomp_x_eval(&pre, yValues, out);
  }
};

// Build degree N and degree 2N precomputes, for multiplication
// For f * g = h, f and g are deg N, h is deg 2N
// and we check the equality by eval on some point
class MultCheckPreComp {
  fmpz_t x;

  const BatchPre* const degN;
  const BatchPre* const deg2N;

  PreX* xN = nullptr;
  PreX* x2N = nullptr;

public:

  MultCheckPreComp(const size_t N)
  : degN(new BatchPre(RootManager(N).getRoots(), N))
  , deg2N(new BatchPre(RootManager(N).getRoots2(), 2 * N))
  {
    if (not isPowerOfTwo(N))
      error_exit("MultCheckPreComp only takes powers of 2");
    fmpz_init(x);
  }

  // Set/update evaluation point
  // since it needs periodic refreshes
  void setEvalPoint(const fmpz_t eval_point) {
    fmpz_set(x, eval_point);

    clearPreX();
    xN = new PreX(degN, x);
    x2N = new PreX(deg2N, x);
  }

  void getEvalPoint(fmpz_t eval_point) const {
    fmpz_set(eval_point, x);
  }

  void Eval(const fmpz_t* const yValues, fmpz_t out) const {
    xN->Eval(yValues, out);
  }

  void Eval2(const fmpz_t* const yValues, fmpz_t out) const {
    x2N->Eval(yValues, out);
  }

  void clearPreX() {
    if (xN != nullptr) {
      delete xN;
      delete x2N;
    }
    xN = nullptr;
    x2N = nullptr;
  }

  ~MultCheckPreComp() {
    clearPreX();
    delete degN;
    delete deg2N;
    fmpz_clear(x);
  }
};


#endif