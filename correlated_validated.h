#ifndef CORRELATED_VALIDATED_H
#define CORRELATED_VALIDATED_H

#include "correlated.h"
#include "interp.h"
#include "utils.h"

class ValidateCorrelatedStore : public PrecomputeStore {
  fmpz_t sigma;
  size_t sigma_uses;
  void check_sigma();  // Call before using sigma. Not return to avoid extra copy.
  void new_sigma();

  const size_t true_batch_size;

  // TODO: merge behavior with precomputed stores
  // Pre: for batch validate
  // unval: recieved
  // val: for computation

  std::queue<const DaBit* const> unvalidated_dabit_store;

  // unvalidated
  std::queue<const AltTriple* const> unvalidated_alt_triple_store;
  std::queue<const AltTriple* const> validated_alt_triple_store;

  std::unordered_map<size_t, MultCheckPreComp*> eval_precomp_store;

  // TODO: consider moving this to normal precompute, and just not using it.
  typedef std::tuple <const DaBit* const *, const AltTriple* const *> pairtype;
  std::unordered_map<std::string, pairtype> unvalidated_pairs;

  // TODO: some way to bulk generate alt triples (e.g. OT).
  const AltTriple* get_validated_alt_triple();

// Batch size must be power of two. NextPowerOfTwo is not inclusive, so -1 to make it so.
public:
  ValidateCorrelatedStore(const int serverfd, const int server_num,
      OT_Wrapper* const ot0, OT_Wrapper* const ot1,
      const size_t batch_size,
      const bool lazy = false)
  : PrecomputeStore(serverfd, server_num, ot0, ot1, NextPowerOfTwo(batch_size-1), lazy)
  , true_batch_size(NextPowerOfTwo(batch_size-1))
  {
    fmpz_init(sigma);

    eval_precomp_store[true_batch_size] = new MultCheckPreComp(true_batch_size);

    new_sigma();
  };

  ~ValidateCorrelatedStore();

  // [z] = x_this * x_other
  // uses a (un)validated alt triple
  int multiplyAltShares(const size_t N, const fmpz_t* const x, fmpz_t* const z,
                        const bool* const use_validated);

  // Add unvalidated correlated to queue
  void addUnvalidated(const DaBit* const dabit, const AltTriple* const trip);
  // Queue up paired unvalidated.
  // Done this way in case data is out of order, had sync issues just using basic add
  // accumulate based on pk, so that it's done with the same "owner"
  void queueUnvalidated(const DaBit* const * dabits, const AltTriple* const * trips,
                        const std::string pk);
  // add up all unvalidated corresponding to pk
  void processUnvalidated(const std::string pk, const size_t n);

  void checkDaBits(const size_t n = 0);
  // void checkTriples(const size_t n = 0);

  int batchValidate(const size_t N);
  int batchValidate() {
    return batchValidate(true_batch_size);
  };
  // TODO: batch validate max possible? NextPowerOfTwo(store size) / 2?

  MultCheckPreComp* getPrecomp(const size_t N);

  void printSizes();

  size_t numvalidated() {
    return dabit_store.size();
  };
};



#endif
