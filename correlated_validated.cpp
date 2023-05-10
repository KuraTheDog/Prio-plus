#include "correlated_validated.h"

#include <sys/wait.h>
#include <iostream>

#include "net_share.h"


const AltTriple* ValidateCorrelatedStore::get_validated_alt_triple() {
  if (validated_alt_triple_store.size() < 1) {
    std::cout << "Doing lazy gen of validated alt triple" << std::endl;
    // TODO: currently lazy gen. Do smart gen. Can do inefficinetly with normal triple.
    AltTriple* t = new AltTriple();
    if (server_num == 0) {
      AltTriple* t_other = new AltTriple();
      NewAltTriples(t, t_other);
      send_AltTriple(serverfd, t_other);
      delete t_other;
    } else {
      recv_AltTriple(serverfd, t);
    }

    return t;
  }
  else {
    const AltTriple* const trip = validated_alt_triple_store.front();
    validated_alt_triple_store.pop();
    return trip;
  }
}

// [z] = x * x'
/*
Trip: (ab, c), (ab', c'), with ab * ab' = c + c'
diff = x - ab, diff' = x' - ab'
swap diff
z_0 = (diff' * trip->AB) + trip->C
z_1 = (diff * x_1) + trip->C
*/
int ValidateCorrelatedStore::multiplyAltShares(
    const size_t N, const fmpz_t* const x, fmpz_t* const z,
    const bool* const use_validated) {
  int sent_bytes = 0;
  fmpz_t* diff; new_fmpz_array(&diff, N);

  // only used by server 1, save trip->AB value
  fmpz_t* a; new_fmpz_array(&a, N);

  for (unsigned int i = 0; i < N; i++) {
    const AltTriple* trip;
    if (use_validated[i]) {
      trip = get_validated_alt_triple();
    } else {
      trip = unvalidated_alt_triple_store.front();
      unvalidated_alt_triple_store.pop();
    }

    fmpz_mod_sub(diff[i], x[i], trip->AB, mod_ctx);
    fmpz_set(z[i], trip->C);
    if (server_num == 0)
      fmpz_set(a[i], trip->AB);
    else
      fmpz_set(a[i], x[i]);

    delete trip;
  }

  // TODO: Swap auto-sums, we want raw swap with threading.
  // sent_bytes += swap_fmpz_batch(serverfd, diff, N);
  sent_bytes += send_fmpz_batch(serverfd, diff, N);
  fmpz_t* diff_other; new_fmpz_array(&diff_other, N);
  recv_fmpz_batch(serverfd, diff_other, N);

  for (unsigned int i = 0; i < N; i++) {
    fmpz_mod_addmul(z[i], diff_other[i], a[i], mod_ctx);
  }
  clear_fmpz_array(diff, N);
  clear_fmpz_array(diff_other, N);
  clear_fmpz_array(a, N);

  return sent_bytes;
}

void ValidateCorrelatedStore::addUnvalidated(
    const DaBit* const dabit, const AltTriple* const trip) {
  unvalidated_dabit_store.push(dabit);
  unvalidated_alt_triple_store.push(trip);
}

void ValidateCorrelatedStore::queueUnvalidated(
    const DaBit* const * dabits, const AltTriple* const * trips,
    const std::string pk) {
  unvalidated_pairs[pk] = {dabits, trips};
}

void ValidateCorrelatedStore::processUnvalidated(const std::string pk, const size_t n) {
  pairtype lists = unvalidated_pairs[pk];
  unvalidated_pairs.erase(pk);
  for (unsigned int i = 0; i < n; i++) {
    addUnvalidated(std::get<0>(lists)[i], std::get<1>(lists)[i]);
  }
  delete[] std::get<0>(lists);
  delete[] std::get<1>(lists);
}

// TODO: Backup process if there is unvalidated.
// Perhaps merge with Preocmpute store, to precompute only if necessary?
// TODO: Use legit triple for validation, for maliciousness

int ValidateCorrelatedStore::batchValidate(const size_t N) {
  int sent_bytes = 0;
  // isPowerOfTwo from utils.h
  if (not isPowerOfTwo(N)) {
    std::cout << N << " should be a power of two" << std::endl;
    return sent_bytes;
  }
  if (unvalidated_dabit_store.size() < N) {
    std::cout << "Not enough unvalidated dabits: " << unvalidated_dabit_store.size();
    std::cout << " of " << N << std::endl;
    return sent_bytes;
  }
  if (unvalidated_dabit_store.size() < N) {
    std::cout << "Not enough unvalidated alt triples: " << unvalidated_alt_triple_store.size();
    std::cout << " of " << N << std::endl;
    return sent_bytes;
  }

  check_sigma();
  sigma_uses += 1;

  MultCheckPreComp* chk = getPrecomp(N);

  const DaBit* candidates[N];

  // Setup
  fmpz_t* pointsF; new_fmpz_array(&pointsF, N);
  fmpz_t* pointsG; new_fmpz_array(&pointsG, 2*N);
  for (unsigned int i = 0; i < N; i++) {
    const DaBit* const bit = unvalidated_dabit_store.front();
    unvalidated_dabit_store.pop();

    // *2 on server0 only
    fmpz_set_ui(pointsF[i], (1 + (server_num == 0)) * (int) bit->b2);
    fmpz_set_ui(pointsG[2*i], (int) bit->b2);
    fmpz_mod_sub(pointsG[2*i], pointsG[2*i], bit->bp, mod_ctx);

    candidates[i] = bit;
  }

  // F on sigma, and extra points
  fmpz_t sigmaF; fmpz_init(sigmaF);
  chk->Eval(pointsF, sigmaF);
  fmpz_t* const polyF = interpolate_inv(N, pointsF);
  clear_fmpz_array(pointsF, N);
  fmpz_t* paddedF; new_fmpz_array(&paddedF, 2*N);
  copy_fmpz_array(paddedF, polyF, N);
  clear_fmpz_array(polyF, N);
  fmpz_t* const evalsF = interpolate_2N(N, paddedF);
  clear_fmpz_array(paddedF, 2*N);

  // Multiply evalF's to get G points
  // Also multiply sigmaF at the same time, as validated
  fmpz_t* points; new_fmpz_array(&points, N+1);
  for (unsigned int i = 0; i < N; i++) {
    fmpz_set(points[i], evalsF[2*i+1]);
  }
  clear_fmpz_array(evalsF, N);
  fmpz_set(points[N], sigmaF);
  // Only the final one needs to be validated.
  // TODO: look back at this logic.
  bool use_validated[N+1];
  memset(use_validated, 0, sizeof(bool)*N);
  use_validated[N] = true;
  fmpz_t* pointsMult; new_fmpz_array(&pointsMult, N+1);
  sent_bytes += multiplyAltShares(N+1, points, pointsMult, use_validated);
  clear_fmpz_array(points, N+1);
  for (unsigned int i = 0; i < N; i++) {
    fmpz_set(pointsG[2*i+1], pointsMult[i]);
  }
  fmpz_set(sigmaF, pointsMult[N]);
  clear_fmpz_array(pointsMult, N+1);

  // Final eval
  // sigmaF0 * sigmaF1 = sigmaG0 + sigmaG1
  // Already did mult, so just add to same, aka differences sum to 0
  fmpz_t sigmaG; fmpz_init(sigmaG);
  chk->Eval2(pointsG, sigmaG);
  clear_fmpz_array(pointsG, 2*N);
  fmpz_t diff; fmpz_init(diff);
  fmpz_mod_sub(diff, sigmaF, sigmaG, mod_ctx);
  fmpz_clear(sigmaF);
  fmpz_clear(sigmaG);

  // Swap, but single item so doens't need wrapper
  fmpz_t diff_other; fmpz_init(diff_other);
  sent_bytes += send_fmpz(serverfd, diff);
  recv_fmpz(serverfd, diff_other);
  fmpz_mod_add(diff, diff, diff_other, mod_ctx);
  fmpz_clear(diff_other);

  if (fmpz_is_zero(diff)) {
    // std::cout << "batch validate is valid" << std::endl;
    for (unsigned int i = 0; i < N; i++) {
      dabit_store.push(candidates[i]);
    }
  } else {
    std::cout << "batch validate is invalid, diff check is " << fmpz_get_si(diff) << std::endl;
    for (unsigned int i = 0; i < N; i++) {
      delete candidates[i];
    }
    error_exit("TODO: Currently no recovery for bad validation");
  }

  fmpz_clear(diff);

  return sent_bytes;
}

void ValidateCorrelatedStore::checkDaBits(const size_t n) {
  if (dabit_store.size() < n) {
    // BatchValidate currently handles the yelling if not enough.
    // NextPowerOfTwo is not inclusive, so -1 or /2

    // minimal approach: smallest needed, not as efficient
    // const size_t to_make = NextPowerOfTwo(n - dabit_store.size() - 1);

    // largest approach: largest power of 2 < size (since exclusive)
    const size_t to_make = NextPowerOfTwo(unvalidated_dabit_store.size())/2;
    batchValidate(to_make);
  }
}

void ValidateCorrelatedStore::check_sigma() {
  if (sigma_uses > EVAL_REUSE_THRESHOLD) {
    new_sigma();
  }
}

void ValidateCorrelatedStore::new_sigma() {
  if (server_num == 0) {
    fmpz_randm(sigma, seed, Int_Modulus);
    send_fmpz(serverfd, sigma);
  } else {
    recv_fmpz(serverfd, sigma);
  }
  sigma_uses = 0;
  for (const auto& pair : eval_precomp_store)
      pair.second -> setEvalPoint(sigma);
}

MultCheckPreComp* ValidateCorrelatedStore::getPrecomp(const size_t N) {
  MultCheckPreComp* pre;
  if (eval_precomp_store.find(N) == eval_precomp_store.end()) {
    pre = new MultCheckPreComp(N);
    pre->setEvalPoint(sigma);
    eval_precomp_store[N] = pre;
  } else {
    pre = eval_precomp_store[N];
  }
  return pre;
}

void ValidateCorrelatedStore::printSizes() {
  std::cout << "Current store sizes:" << std::endl;
  std::cout << " Dabits: " << dabit_store.size() << std::endl;
  std::cout << "  Unvalidated: " << unvalidated_dabit_store.size() << std::endl;
  std::cout << " Alt Triples: " << validated_alt_triple_store.size() << std::endl;
  std::cout << "  Unvalidated: " << unvalidated_alt_triple_store.size() << std::endl;
}

ValidateCorrelatedStore::~ValidateCorrelatedStore() {
  fmpz_clear(sigma);
  for (const auto& pair : eval_precomp_store)
    delete pair.second;

  while (!unvalidated_dabit_store.empty()) {
    const DaBit* const bit = unvalidated_dabit_store.front();
    unvalidated_dabit_store.pop();
    delete bit;
  }
  while (!validated_alt_triple_store.empty()) {
    const AltTriple* const trip = validated_alt_triple_store.front();
    validated_alt_triple_store.pop();
    delete trip;
  }
  while (!unvalidated_alt_triple_store.empty()) {
    const AltTriple* const trip = unvalidated_alt_triple_store.front();
    unvalidated_alt_triple_store.pop();
    delete trip;
  }
}
