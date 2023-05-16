#include <iostream>

#include "constants.h"
#include "fmpz_utils.h"
#include "interp.h"

extern "C" {
    #include "poly/fft.h"
}

std::unordered_map<size_t, fmpz_t*> RootManager::roots;
std::unordered_map<size_t, fmpz_t*> RootManager::roots2;
std::unordered_map<size_t, fmpz_t*> RootManager::invroots;
std::unordered_map<size_t, fmpz_t*> RootManager::invroots2;

void RootManager::addRoots(const size_t N){
  if (roots.count(N) == 1)
    return;

  fmpz_t* r; new_fmpz_array(&r, N);
  fmpz_t* inv; new_fmpz_array(&inv, N);
  fmpz_t* r2; new_fmpz_array(&r2, 2 * N);
  fmpz_t* r2inv; new_fmpz_array(&r2inv, 2 * N);

  const int step_size = (1 << twoOrder) / N;
  fmpz_t g_; fmpz_init(g_);          // Generator order N
  fmpz_t ghalf_; fmpz_init(ghalf_);  // g_^(step/2), for 2N roots

  /*
    N = 2^k, so stepsize = 2^(Ord - k).
    g_ = gen^stepsize.
    So g_^N = gen^(2^ord) = 1, by fermat little.
  */
  fmpz_powm_ui(g_, Int_Gen, step_size, Int_Modulus);
  fmpz_powm_ui(ghalf_, Int_Gen, step_size / 2, Int_Modulus);

  fmpz_set_ui(r[0], 1);
  fmpz_set_ui(inv[0], 1);
  fmpz_set_ui(r2[0], 1);
  fmpz_set_ui(r2inv[0], 1);

  for (unsigned int i = 1; i < N; i++) {
    fmpz_mod_mul(r[i], r[i-1], g_, mod_ctx);
    fmpz_set(inv[N-i], r[i]);
  }
  for (unsigned int i = 1; i < 2 * N; i++) {
    fmpz_mod_mul(r2[i], r2[i-1], ghalf_, mod_ctx);
    fmpz_set(r2inv[2*N-i], r2[i]);
  }

  fmpz_clear(g_);
  fmpz_clear(ghalf_);

  roots[N] = r;
  invroots[N] = inv;
  roots2[N] = r2;
  invroots2[N] = r2inv;
}

fmpz_t* interpolate_N(const size_t N, const fmpz_t* const points) {
  return fft_interpolate(Int_Modulus, N, RootManager(N).getRoots(), points, false);
}

fmpz_t* interpolate_N_inv(const size_t N, const fmpz_t* const points) {
  return fft_interpolate(Int_Modulus, N, RootManager(N).getRootsInv(), points, true);
}

fmpz_t* interpolate_2N(const size_t N, const fmpz_t* const points) {
  return fft_interpolate(Int_Modulus, 2 * N, RootManager(N).getRoots2(), points, false);
}

fmpz_t* interpolate_2N_inv(const size_t N, const fmpz_t* const points) {
  return fft_interpolate(Int_Modulus, 2 * N, RootManager(N).getRoots2Inv(), points, true);
}

void MultEvalManager::check_eval_point(const size_t n) {
  if (eval_point_uses > EVAL_REUSE_THRESHOLD) {
    new_eval_point();
  }
  eval_point_uses += n;
}

void MultEvalManager::new_eval_point() {
  if (server_num == 0) {
    fmpz_randm(eval_point, seed, Int_Modulus);
    send_fmpz(serverfd, eval_point);
  } else {
    recv_fmpz(serverfd, eval_point);
  }
  eval_point_uses = 0;
  for (const auto& pair : store)
    pair.second -> setEvalPoint(eval_point);
}

MultCheckPreComp* MultEvalManager::get_Precomp(const size_t N) {
  if (store.find(N) == store.end()) {
    store[N] = new MultCheckPreComp(N);
    store[N]->setEvalPoint(eval_point);
  }
  return store[N];
}
