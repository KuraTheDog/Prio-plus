#ifndef CORRELATED_VALIDATED_H
#define CORRELATED_VALIDATED_H

#include "correlated.h"
#include "interp.h"
#include "utils.h"

class ValidateCorrelatedStore : public DaBitStore {
  fmpz_t sigma;
  size_t sigma_uses;
  void check_sigma();  // Call before using sigma. Not return to avoid extra copy.
  void new_sigma();

  std::queue<const DaBit* const> unvalidated_dabit_store;

  // unvalidated
  std::queue<const AltTriple* const> unvalidated_alt_triple_store;
  std::queue<const AltTriple* const> validated_alt_triple_store;

  MultCheckPreComp* chk;

  // TODO: some way to bulk generate alt triples (e.g. OT).
  const AltTriple* get_validated_alt_triple();

// Batch size must be power of two. NextPowerOfTwo is not inclusive, so -1 to make it so.
public:
  ValidateCorrelatedStore(const int serverfd, const int server_num,
      const size_t batch_size,
      const bool lazy = false, const bool do_fork = true)
  : DaBitStore(serverfd, server_num, NextPowerOfTwo(batch_size-1), lazy, do_fork)
  , chk(new MultCheckPreComp(NextPowerOfTwo(batch_size-1)))
  {
    fmpz_init(sigma);

    new_sigma();
  };

  ~ValidateCorrelatedStore();

  fmpz_t* multiplyAltShares(const size_t N, const fmpz_t* const x, 
                            const bool* const use_validated);

  void addUnvalidated(const DaBit* const dabit, const AltTriple* const trip);

  void checkDaBits(const size_t n = 0);

  void batchValidate(const size_t N);
  void batchValidate() {
    batchValidate(batch_size);
  };
  // TODO: batch validate max possible? NextPowerOfTwo(store size) / 2? 

  size_t numvalidated() {
    return dabit_store.size();
  };
};



#endif
