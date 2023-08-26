#include "correlated_validated.h"

#include <sys/wait.h>
#include <iostream>

#include "net_share.h"

const AltTriple* const * const ValidateCorrelatedStore::gen_AltTriple(const size_t N) {
  // TODO: smart gen. Can build out of normal beaver triples.
  std::cout << "TODO: Currently only lazy alt triple gen\n";
  return gen_AltTriple_lazy(N);
}

const AltTriple* const * const ValidateCorrelatedStore::gen_AltTriple_lazy(const size_t N) {
  AltTriple** const t = new AltTriple*[N];
  for (unsigned int i = 0; i < N; i++)
    t[i] = new AltTriple();

  if (server_num == 0) {
    AltTriple** const t_other = new AltTriple*[N];
    for (unsigned int i = 0; i < N; i++) {
      t_other[i] = new AltTriple();
      makeLocalAltTriple(t[i], t_other[i]);
    }
    send_AltTriple_batch(serverfd, t_other, N);
    for (unsigned int i = 0; i < N; i++)
      delete t_other[i];
    delete[] t_other;
  } else {
    recv_AltTriple_batch(serverfd, t, N);
  }
  return t;
}

const AltTriple* ValidateCorrelatedStore::get_AltTriple(const bool validated) {
  check_AltTriple(1, validated);

  if (!validated and unvalidated_alt_triple_store.size() > 0) {
    const AltTriple* const t = unvalidated_alt_triple_store.front();
    unvalidated_alt_triple_store.pop();
    return t;
  }

  const AltTriple* const t = alt_triple_store.front();
  alt_triple_store.pop();
  return t;
}

void ValidateCorrelatedStore::check_AltTriple(const size_t n, const bool validated) {
  size_t num = alt_triple_store.size();
  if (!validated)
    num += unvalidated_alt_triple_store.size();
  if (num < n) {
    add_AltTriples(n - num);
  }
}

void ValidateCorrelatedStore::check_AltTriple(const size_t n, const bool* const validated) {
  // std::cout << "Calling mixed check_AltTriples(" << n << ")\n";
  size_t val_count = 0;
  for (unsigned int i = 0; i < n; i++)
    val_count += (size_t) validated[i];

  // If not enough unval (need n - val), then more val
  if (unvalidated_alt_triple_store.size() < n - val_count) {
    val_count += (n - val_count) - unvalidated_alt_triple_store.size();
  }
  if (alt_triple_store.size() < val_count) {
    add_AltTriples(val_count - alt_triple_store.size());
  }
}

void ValidateCorrelatedStore::add_AltTriples(const size_t n) {
  auto start = clock_start();
  const size_t num_to_make = (n > alt_triple_batch_size ? n : alt_triple_batch_size);

  std::cout << "adding " << num_to_make << " AltTriples" << std::endl;
  const AltTriple* const * t;

  if (lazy)
    t = gen_AltTriple_lazy(num_to_make);
  else
    t = gen_AltTriple(num_to_make);

  for (unsigned int i = 0; i < num_to_make; i++)
    alt_triple_store.push(t[i]);
  delete[] t;

  std::cout << "add_AltTriples timing : " << sec_from(start) << std::endl;
}

// [z] = x * x'
/*
Trip: (ab, c), (ab', c'), with ab * ab' = c + c'
diff = x - ab, diff' = x' - ab'
swap diff
z_0 = (diff' * trip->AB) + trip->C
z_1 = (diff * x_1) + trip->C
*/
int64_t ValidateCorrelatedStore::multiply_AltShares(
    const size_t N, const fmpz_t* const x, fmpz_t* const z,
    const bool* const validated) {
  int64_t sent_bytes = 0;
  fmpz_t* diff; new_fmpz_array(&diff, N);

  // only used by server 1, save trip->AB value
  fmpz_t* a; new_fmpz_array(&a, N);

  check_AltTriple(N, validated);

  for (unsigned int i = 0; i < N; i++) {
    const AltTriple* trip = get_AltTriple(validated[i]);

    fmpz_mod_sub(diff[i], x[i], trip->AB, mod_ctx);
    fmpz_set(z[i], trip->C);
    if (server_num == 0)
      fmpz_set(a[i], trip->AB);
    else
      fmpz_set(a[i], x[i]);

    delete trip;
  }

  fmpz_t* diff_other; new_fmpz_array(&diff_other, N);
  sent_bytes += swap_fmpz_batch(serverfd, diff, diff_other, N);

  for (unsigned int i = 0; i < N; i++) {
    fmpz_mod_addmul(z[i], diff_other[i], a[i], mod_ctx);
  }
  clear_fmpz_array(diff, N);
  clear_fmpz_array(diff_other, N);
  clear_fmpz_array(a, N);

  return sent_bytes;
}

void ValidateCorrelatedStore::add_Unvalidated(
    const DaBit* const dabit, const AltTriple* const trip) {
  unvalidated_dabit_store.push(dabit);
  unvalidated_alt_triple_store.push(trip);
}

void ValidateCorrelatedStore::queue_Unvalidated(
    const DaBit* const * dabits, const AltTriple* const * trips,
    const std::string tag, const size_t n) {
  unvalidated_pairs[tag] = {n, dabits, trips};
}

void ValidateCorrelatedStore::process_Unvalidated(const std::string tag) {
  pairtype lists = unvalidated_pairs[tag];
  unvalidated_pairs.erase(tag);
  const size_t n = std::get<0>(lists);
  for (unsigned int i = 0; i < n; i++) {
    add_Unvalidated(std::get<1>(lists)[i], std::get<2>(lists)[i]);
  }
  delete[] std::get<1>(lists);
  delete[] std::get<2>(lists);
}

void ValidateCorrelatedStore::batch_Validate_param(
    const size_t target, size_t& N, size_t& num_val) {
  const size_t target_2 = NextPowerOfTwo(target - 1);
  const size_t num_unval = unvalidated_dabit_store.size();
  if (target_2 <= num_unval) {
    // Make as many as possible within unval
    // Largest power of 2 <= num_unval
    N = NextPowerOfTwo(num_unval) / 2;
    num_val = N;
  } else {
    // Use up all unvalidated, and also use dummies
    // Smallest power of 2 >= num_unval
    N = NextPowerOfTwo(num_unval - 1);
    num_val = num_unval;
  }
}

void ValidateCorrelatedStore::batch_Validate_setup(
    const size_t N, const size_t num_val, fmpz_t* const points, fmpz_t* const pointsMult,
    bool* const use_validated, const DaBit** const candidates) {
  mult_eval_manager.check_eval_point(2);
  // Should always be synced, but just in case. Maybe complain instead?
  check_AltTriple(N, false);
  // Also have a normal validated one, i.e. from outside
  check_AltTriple(1, true);

  // Setup polynomial points
  // Note that new_fmpz_array zero initializes, which give "dummy" dabits.
  fmpz_t* pointsF; new_fmpz_array(&pointsF, N);
  for (unsigned int i = 0; i < num_val; i++) {
    const DaBit* const bit = unvalidated_dabit_store.front();
    unvalidated_dabit_store.pop();

    // *2 on server0 only
    fmpz_set_ui(pointsF[i], (1 + (server_num == 0)) * (int) bit->b2);

    candidates[i] = bit;
  }

  // F on sigma, and extra points
  fmpz_t sigmaF; fmpz_init(sigmaF);
  mult_eval_manager.get_Precomp(N)->Eval(pointsF, sigmaF);
  fmpz_t* const polyF = interpolate_N_inv(N, pointsF);
  clear_fmpz_array(pointsF, N);
  fmpz_t* paddedF; new_fmpz_array(&paddedF, 2*N);
  copy_fmpz_array(paddedF, polyF, N);
  clear_fmpz_array(polyF, N);
  fmpz_t* const evalsF = interpolate_2N(N, paddedF);
  clear_fmpz_array(paddedF, 2*N);

  // Multiply evalF's to get G points
  // Also multiply sigmaF at the same time, as validated
  for (unsigned int i = 0; i < N; i++)
    fmpz_set(points[i], evalsF[2*i+1]);
  clear_fmpz_array(evalsF, N);
  fmpz_set(points[N], sigmaF);
  // Only the final one needs to be validated.
  memset(use_validated, 0, sizeof(bool)*N);
  use_validated[N] = true;
}

void ValidateCorrelatedStore::batch_Validate_process(
    const size_t N, const size_t num_val, const fmpz_t* const pointsMult,
    const DaBit** const candidates, fmpz_t diff
  ) {
  // Step 2: Process and evaluate
  fmpz_t* pointsG; new_fmpz_array(&pointsG, 2*N);
  for (unsigned int i = 0; i < num_val; i++) {
    fmpz_set_ui(pointsG[2*i], (int) candidates[i]->b2);
    fmpz_mod_sub(pointsG[2*i], pointsG[2*i], candidates[i]->bp, mod_ctx);
  }
  for (unsigned int i = 0; i < N; i++) {
    fmpz_set(pointsG[2*i+1], pointsMult[i]);
  }
  fmpz_t sigmaF; fmpz_init_set(sigmaF, pointsMult[N]);

  // Final eval
  // sigmaF0 * sigmaF1 = sigmaG0 + sigmaG1
  // Already did mult, so just add to same, aka differences sum to 0
  fmpz_t sigmaG; fmpz_init(sigmaG);
  mult_eval_manager.get_Precomp(N)->Eval2(pointsG, sigmaG);
  clear_fmpz_array(pointsG, 2*N);
  fmpz_mod_sub(diff, sigmaF, sigmaG, mod_ctx);
  fmpz_clear(sigmaF);
  fmpz_clear(sigmaG);
}

void ValidateCorrelatedStore::batch_Validate_finish(
    const size_t num_val, fmpz_t diff, const DaBit** const candidates) {
  if (fmpz_is_zero(diff)) {
    // std::cout << "batch validate is valid" << std::endl;
    // Only push those corresponding to original unvalidated
    for (unsigned int i = 0; i < num_val; i++) {
      validated_dabit_store.push(candidates[i]);
    }
  } else {
    std::cout << "batch validate is invalid, diff check is " << fmpz_get_si(diff) << std::endl;
    for (unsigned int i = 0; i < num_val; i++) {
      delete candidates[i];
    }
    // Assuming called through check_DaBits, it will then proceed to use precomputed daBits instead.
  }
}

int64_t ValidateCorrelatedStore::batch_Validate(const size_t target) {
  int64_t sent_bytes = 0;

  // Step 0: Setup

  // Param setup
  size_t N;        // num to make, power of 2, including dummies
  size_t num_val;  // how many to validate
  batch_Validate_param(target, N, num_val);

  // std::cout << "batch Val post check sizes:\n";
  // print_Sizes();

  const DaBit** const candidates = new const DaBit*[num_val];
  fmpz_t* points; new_fmpz_array(&points, N+1);
  fmpz_t* pointsMult; new_fmpz_array(&pointsMult, N+1);
  bool use_validated[N+1];
  batch_Validate_setup(N, num_val, points, pointsMult, use_validated, candidates);

  // Step 1: Both multiplies
  sent_bytes += multiply_AltShares(N+1, points, pointsMult, use_validated);

  clear_fmpz_array(points, N+1);
  fmpz_t diff; fmpz_init(diff);
  batch_Validate_process(N, num_val, pointsMult, candidates, diff);
  clear_fmpz_array(pointsMult, N+1);

  // Step 3: Reveal differences
  sent_bytes += reveal_fmpz(serverfd, diff);

  batch_Validate_finish(num_val, diff, candidates);

  delete[] candidates;
  fmpz_clear(diff);

  return sent_bytes;
}

int64_t ValidateCorrelatedStore::check_DaBits(const size_t n) {
  int64_t sent_bytes = 0;
  // std::cout << "validated check_DaBits(" << n << ") vs " << num_validated_dabits() << std::endl;
  // If not enough validated, validate everything
  if (num_validated_dabits() < n) {
    // std::cout << "batch validating: " << n - num_validated_dabits() << std::endl;
    sent_bytes += batch_Validate(n - num_validated_dabits());
  }
  // If still not enough validated, ensure enough precomputed
  if (num_validated_dabits() < n) {
    int num_to_make = n - num_validated_dabits();
    // std::cout << "precompute checking: " << n - num_validated_dabits() << std::endl;
    // std::cout << "check_da(" << n << "), not enough validated. to_make = " << num_to_make << "\n";
    sent_bytes += PrecomputeStore::check_DaBits(num_to_make);
    check_AltTriple(num_to_make, true);
  }
  // printSizes();
  return sent_bytes;
}

const DaBit* const ValidateCorrelatedStore::get_DaBit() {
  check_DaBits(1);
  // Default to precomputed if not enough valids
  // std::cout << "getting a " << (num_validated_dabits() > 1 ? "validated" : "precomputed");
  // std::cout << " dabit" << std::endl;
  std::queue<const DaBit*>& q = (
    (num_validated_dabits() > 0) ? validated_dabit_store : dabit_store);
  const DaBit* const ans = q.front();
  q.pop();
  return ans;
}

void ValidateCorrelatedStore::print_Sizes() const {
  std::cout << "Current store sizes:\n";
  std::cout << " Dabits: \n";
  std::cout << "  Precomputed: " << dabit_store.size() << "\n";
  std::cout << "  Unvalidated: " << unvalidated_dabit_store.size() << "\n";
  std::cout << "  Validated:   " << validated_dabit_store.size() << "\n";
  std::cout << " Alt Triples: \n";
  std::cout << "  Unvalidated: " << unvalidated_alt_triple_store.size() << "\n";
  std::cout << "  Validated:   " << alt_triple_store.size() << "\n";
  std::cout << " unvalidated_pairs: " << unvalidated_pairs.size() << std::endl;
}

void ValidateCorrelatedStore::maybe_Update() {
  PrecomputeStore::maybe_Update();
  check_AltTriple(batch_size, true);
}
