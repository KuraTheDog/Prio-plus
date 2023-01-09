#include "correlated_validated.h"

#include <sys/wait.h>
#include <iostream>

#include "net_share.h"


const AltTriple* ValidateCorrelatedStore::get_validated_alt_triple() {
  if (validated_alt_triple_store.size() < 1) {
    std::cout << "Doing lazy gen of validated alt triple" << std::endl;
    // TODO: currently lazy gen. Do smart gen
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

// [z] = x_this * x_other
/* 
diff_i = x-i - trip_i->AB
swap diff
z_0 = (diff_1 * trip->AB) + trip->C
z_1 = (diff_0 * x_1) + trip->C
*/
fmpz_t* ValidateCorrelatedStore::multiplyAltShares(
    const size_t N, const fmpz_t* const x, const bool* const use_validated) {
  fmpz_t* z; new_fmpz_array(&z, N);
  fmpz_t* diff; new_fmpz_array(&diff, N);

  // only used by server 1, save trip->AB value 
  fmpz_t* a; if (server_num == 0) new_fmpz_array(&a, N);

  for (unsigned int i = 0; i < N; i++) {
    const AltTriple* trip;
    if (use_validated[i]) {
      trip = get_validated_alt_triple();
    } else {
      trip = unvalidated_alt_triple_store.front();
      unvalidated_alt_triple_store.pop();
    }

    fmpz_sub(diff[i], x[i], trip->AB);
    fmpz_mod(diff[i], diff[i], Int_Modulus);
    fmpz_set(z[i], trip->C);
    if (server_num == 0) fmpz_set(a[i], trip->AB);

    delete trip;
  }

  pid_t pid = 0;
  int status = 0;
  if (do_fork) pid = fork();
  if (pid == 0) {
    send_fmpz_batch(serverfd, diff, N);
    if (do_fork) exit(EXIT_SUCCESS);
  }
  fmpz_t* diff_other; new_fmpz_array(&diff_other, N);
  recv_fmpz_batch(serverfd, diff_other, N);

  for (unsigned int i = 0; i < N; i++) {
    if (server_num == 0) {
      fmpz_addmul(z[i], diff_other[i], a[i]);
    } else {
      fmpz_addmul(z[i], diff_other[i], x[i]);
    }
    fmpz_mod(z[i], z[i], Int_Modulus);
  }
  clear_fmpz_array(diff, N);
  clear_fmpz_array(diff_other, N);
  if (server_num == 0) clear_fmpz_array(a, N);

  if (do_fork) waitpid(pid, &status, 0);

  return z;
}

void ValidateCorrelatedStore::addUnvalidated(
    const DaBit* const dabit, const AltTriple* const trip) {
  unvalidated_dabit_store.push(dabit);
  unvalidated_alt_triple_store.push(trip);
}

void ValidateCorrelatedStore::batchValidate(const size_t N) {
  // isPowerOfTwo from utils.h
  if (not isPowerOfTwo(N)) {
    std::cout << N << " should be a power of two" << std::endl;
    return;
  }
  if (unvalidated_dabit_store.size() < N) {
    std::cout << "Not enough unvalidated dabits: " << unvalidated_dabit_store.size();
    std::cout << " of " << N << std::endl;
    return;
  }
  if (unvalidated_dabit_store.size() < N) {
    std::cout << "Not enough unvalidated alt triples: " << unvalidated_alt_triple_store.size();
    std::cout << " of " << N << std::endl;
    return;
  }

  check_sigma();
  sigma_uses += 1;

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
    fmpz_sub(pointsG[2*i], pointsG[2*i], bit->bp);
    fmpz_mod(pointsG[2*i], pointsG[2*i], Int_Modulus);

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
  bool use_validated[N+1];
  memset(use_validated, 0, sizeof(bool)*N);
  use_validated[N] = true;
  fmpz_t* pointsMult = multiplyAltShares(N+1, points, use_validated);
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
  fmpz_sub(diff, sigmaF, sigmaG);
  fmpz_mod(diff, diff, Int_Modulus);
  fmpz_t diff_other; fmpz_init(diff_other);

  pid_t pid = 0;
  int status = 0;
  if (do_fork) pid = fork();
  if (pid == 0) {
    send_fmpz(serverfd, diff);
    if (do_fork) exit(EXIT_SUCCESS);
  }
  recv_fmpz(serverfd, diff_other);

  fmpz_add(diff, diff, diff_other);
  fmpz_mod(diff, diff, Int_Modulus);

  fmpz_clear(sigmaF);
  fmpz_clear(sigmaG);
  fmpz_clear(diff_other);

  if (fmpz_is_zero(diff)) {
    for (unsigned int i = 0; i < N; i++) {
      dabit_store.push(candidates[i]);
    }
  } else {
    for (unsigned int i = 0; i < N; i++) {
      delete candidates[i];
    }
  }
  
  fmpz_clear(diff);

  if (do_fork) waitpid(pid, &status, 0);
}

void ValidateCorrelatedStore::checkDaBits(const size_t n) {
  if (dabit_store.size() < n) {
    // Can fail if not enough unvalidated, in which case batchValidate complains
    batchValidate(NextPowerOfTwo(n - dabit_store.size()));
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
  chk->setEvalPoint(sigma);
}

ValidateCorrelatedStore::~ValidateCorrelatedStore() {
  fmpz_clear(sigma);
  delete chk;

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
