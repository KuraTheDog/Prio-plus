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

  std::queue<const AltTriple*>& q = validated ? alt_triple_store : unvalidated_alt_triple_store;

  const AltTriple* const t = q.front();
  q.pop();
  return t;
}

void ValidateCorrelatedStore::check_AltTriple(const size_t n, const bool validated) {
  std::queue<const AltTriple*>& q = validated ? alt_triple_store : unvalidated_alt_triple_store;
  // std::cout << "checking " << n << (validated?" ":" un");
  // std::cout << "validated AltTriples, have: " << q.size() << std::endl;
  if (q.size() < n) {
    // std::cout << "Calling check_AltTriples(" << n << "), only have " << q.size() << "\n";
    add_AltTriples(n - q.size(), validated);
  }
}

void ValidateCorrelatedStore::check_AltTriple(const size_t n, const bool* const validated) {
  // std::cout << "Calling mixed check_AltTriples(" << n << ")\n";
  size_t val_count = 0;
  for (unsigned int i = 0; i < n; i++)
    val_count += (size_t) validated[i];
  check_AltTriple(val_count, true);
  check_AltTriple(n - val_count, false);
}

void ValidateCorrelatedStore::add_AltTriples(const size_t n, const bool validated) {
  auto start = clock_start();
  // std::cout << "Calling add_AltTriples(" << n << ")\n";
  const size_t num_to_make = (n > alt_triple_batch_size ? n : alt_triple_batch_size);

  std::cout << "adding " << num_to_make << (validated?" ":" un");
  std::cout << "validated AltTriples" << std::endl;
  const AltTriple* const * t;
  std::queue<const AltTriple*>& q = validated ? alt_triple_store : unvalidated_alt_triple_store;

  if (lazy)
    t = gen_AltTriple_lazy(num_to_make);
  else
    t = gen_AltTriple(num_to_make);

  for (unsigned int i = 0; i < num_to_make; i++)
    q.push(t[i]);
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

int64_t ValidateCorrelatedStore::batch_Validate(const size_t target) {
  int64_t sent_bytes = 0;

  // Step 0: Setup

  // Param setup
  const size_t target_2 = NextPowerOfTwo(target - 1);
  const size_t num_unval = unvalidated_dabit_store.size();
  size_t N;        // num to make, power of 2, including dummies
  size_t num_val;  // how many to validate
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

  // Should always be synced, but just in case. Maybe complain instead?
  check_AltTriple(N, false);

  // More variable setup
  mult_eval_manager.check_eval_point(2);
  MultCheckPreComp* chk = mult_eval_manager.get_Precomp(N);

  const DaBit* candidates[num_val];

  // Setup polynomial points
  // Note that new_fmpz_array zero initializes, which give "dummy" dabits.
  fmpz_t* pointsF; new_fmpz_array(&pointsF, N);
  fmpz_t* pointsG; new_fmpz_array(&pointsG, 2*N);
  for (unsigned int i = 0; i < num_val; i++) {
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
  fmpz_t* const polyF = interpolate_N_inv(N, pointsF);
  clear_fmpz_array(pointsF, N);
  fmpz_t* paddedF; new_fmpz_array(&paddedF, 2*N);
  copy_fmpz_array(paddedF, polyF, N);
  clear_fmpz_array(polyF, N);
  fmpz_t* const evalsF = interpolate_2N(N, paddedF);
  clear_fmpz_array(paddedF, 2*N);

  // Multiply evalF's to get G points
  // Also multiply sigmaF at the same time, as validated
  fmpz_t* points; new_fmpz_array(&points, N+1);
  for (unsigned int i = 0; i < N; i++)
    fmpz_set(points[i], evalsF[2*i+1]);
  clear_fmpz_array(evalsF, N);
  fmpz_set(points[N], sigmaF);
  // Only the final one needs to be validated.
  bool use_validated[N+1];
  memset(use_validated, 0, sizeof(bool)*N);
  use_validated[N] = true;
  fmpz_t* pointsMult; new_fmpz_array(&pointsMult, N+1);

  // Step 1: Both multiplies
  sent_bytes += multiply_AltShares(N+1, points, pointsMult, use_validated);

  // Step 2: Process and evaluate
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

  // Step 3: Reveal differences
  fmpz_t diff_other; fmpz_init(diff_other);
  sent_bytes += send_fmpz(serverfd, diff);
  recv_fmpz(serverfd, diff_other);

  // Step 4: Check if differences sum to 0
  fmpz_mod_add(diff, diff, diff_other, mod_ctx);
  fmpz_clear(diff_other);

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
  std::cout << " Unvalidated Alt Triples: " << unvalidated_alt_triple_store.size() << std::endl;
}
