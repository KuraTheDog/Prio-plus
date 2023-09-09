#ifndef CORRELATED_VALIDATED_H
#define CORRELATED_VALIDATED_H

#include "correlated.h"
#include "interp.h"
#include "utils.h"

class ValidateCorrelatedStore : public PrecomputeStore {
protected:
  const size_t min_batch_size;
  const size_t alt_triple_batch_size;

  // Pre: for batch validate
  // unval: received
  // val: for computation
  /* Inherited
  std::queue<const DaBit*> dabit_store;
  std::queue<const BooleanBeaverTriple*> btriple_store;
  std::queue<const BeaverTriple*> atriple_store;
  */

  // TODO: Add in arith triples. Maybe bool, but those currently unused.
  std::queue<const DaBit*> unvalidated_dabit_store;
  std::queue<const DaBit*> validated_dabit_store;
  // Precomputed (inherited)
  // std::queue<const DaBit*> dabit_store;

  std::queue<const AltTriple*> unvalidated_alt_triple_store;
  // Precomputed
  std::queue<const AltTriple*> alt_triple_store;

  MultEvalManager mult_eval_manager;

  typedef std::tuple <size_t, const DaBit* const *, const AltTriple* const *> pairtype;
  std::unordered_map<std::string, pairtype> unvalidated_pairs;

  // TODO: some way to bulk generate alt triples (e.g. OT).
  // If unvalidated, can draw from validated as backup
  const AltTriple* get_AltTriple(const bool validated);
  const AltTriple* const * const gen_AltTriple(const size_t N);
  const AltTriple* const * const gen_AltTriple_lazy(const size_t N);
  void add_AltTriples(const size_t n);

  const DaBit* const get_DaBit() override;

// Batch size must be power of two. NextPowerOfTwo is not inclusive, so -1 to make it so.
public:
  ValidateCorrelatedStore(const int serverfd, const int server_num,
      OT_Wrapper* const ot0, OT_Wrapper* const ot1,
      const size_t batch_size,
      const int lazy = 0)
  : PrecomputeStore(serverfd, server_num, ot0, ot1, batch_size, lazy)
  , min_batch_size(NextPowerOfTwo(batch_size-1))
  , alt_triple_batch_size(batch_size)
  , mult_eval_manager(server_num, serverfd)
  {
    // Common, so pre-cache
    mult_eval_manager.get_Precomp(min_batch_size);
  };

  ~ValidateCorrelatedStore() {
    clear_queue(unvalidated_dabit_store);
    clear_queue(validated_dabit_store);
    clear_queue(unvalidated_alt_triple_store);
    clear_queue(alt_triple_store);

    for (auto it = unvalidated_pairs.begin(); it != unvalidated_pairs.end(); ++it) {
      pairtype p = it->second;
      for (unsigned int i = 0; i < std::get<0>(p); i++) {
        delete std::get<1>(p)[i];
        delete std::get<2>(p)[i];
      }
      delete[] std::get<1>(p);
      delete[] std::get<2>(p);
    }
    unvalidated_pairs.clear();
  }

  // if unvalidated, checks both stores.
  // If validated, checks only validated store
  void check_AltTriple(const size_t n, const bool validated);
  // Takes in list for convenience. just checks for #valid and #unvalid
  void check_AltTriple(const size_t n, const bool* const validated);

  // [z] = x_this * x_other
  // uses a (un)validated alt triple
  void multiply_AltShares_setup(
      const size_t N, const bool* const validated, const fmpz_t* const x,
      fmpz_t* const z, fmpz_t* const diff, fmpz_t* const a);
  int64_t multiply_AltShares(const size_t N,
      const fmpz_t* const x, fmpz_t* const z, const bool* const validated);

  // Queue up paired unvalidated.
  // Done this way in case data is out of order between the two servers
  // accumulate based on tag, so that shares are aligned
  // Arrays of pointers each length n.
  void queue_Unvalidated(const DaBit* const * dabits, const AltTriple* const * trips,
                         const std::string tag, const size_t n);
  // add in all n unvalidated corresponding to tag, where n is the size of the lists queued
  // For use when generally iterating over (synced) tags.
  void process_Unvalidated(const std::string tag);
  // Add unvalidated correlated to queue
  void add_Unvalidated(const DaBit* const dabit, const AltTriple* const trip);

  int64_t check_DaBits(const size_t n = 0) override;
  // int64_t checkTriples(const size_t n = 0);

  /*
  Tries to make at least N validated dabits.
    N' = NextPowerOf2(N) inclusive
  Depends on how many num_unvalid there are.
  a) N' <= num_unvalid
      Then it can make the largest power of 2 it can using only unvalid.
  b) N' > num_unvalid
      Uses "dummy" zero dabits to pad unvalid dabits.
      Pulls extra alt-triples from precomputed.

  Loosely:
  Takes N unvalidated DaBits and N' unvalidated altTriples
  Also uses 1 precomputed altTrip
  Produces N validated daBits

  2 rounds:
  Round 1: Alt multiplication
    N unvalidated for extension
    Simultaneously 1 validated for on the eval point sigma
  Round 2: Reveal
    Swap computed (fg - h)(sigma)
  */
  // N: num to make, power of 2, including dummies, num_val: how many to validate
  void batch_Validate_param(const size_t target, size_t& N, size_t& num_val);
  void batch_Validate_setup(const size_t N, const size_t num_val, fmpz_t* const pointsMult,
      fmpz_t* const delta, fmpz_t* const a, const DaBit** const candidates);
  void batch_Validate_process(const size_t N, const size_t num_val, fmpz_t* const pointsMult,
      const fmpz_t* const delta_other, const fmpz_t* const a, const DaBit** const candidates, fmpz_t diff);
  void batch_Validate_finish(const size_t num_val, fmpz_t diff, const DaBit** const candidates);
  int64_t batch_Validate(const size_t N);
  int64_t batch_Validate() {
    return batch_Validate(min_batch_size);
  }

  MultCheckPreComp* get_Precomp(const size_t N);

  void print_Sizes() const override;
  void maybe_Update() override; // Precompute if not enough.

  size_t num_validated_dabits() const {
    return validated_dabit_store.size();
  };
};

#endif
