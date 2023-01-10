#ifndef CORRELATED_H
#define CORRELATED_H

/*
Correlated precomputes
synced correlated precomputed values between servers
Maintaines a cache of precomputed correlated randomness.

Due to send buffers potentially filling up,
send is done on a forked child, while parent receives
This can be turned off for testing.

PrecomputeStore:
Made in at least batch_size chunks at a time.
New batches are build either as it runs out, or by calling maybeUpdate
Also includes io objects for OT, which have their own precomputes
Uses for generating daBits via OT.
Also adding boolean shares via boolean triples (currently unused).
OT share conversion only for comparison with USE_OT_B2A.

ValidateCorrelatedStore:
Clients also send unvalidated randomness, then it batch validates them.
Has validated and unvalidated cache. Uses validated for compute, while unvalided fill up.

OTCorrelatedStore:
Pure OT based B2A, via intsum_ot
*/

#include <emp-ot/emp-ot.h>
#include <emp-tool/emp-tool.h>
#include <queue>

#include "ot.h"
#include "share.h"

// Selector
enum StoreType {
  precompute,
  ot,
};

// Global for now. Can have client/server sync it or be a param file, but not necessary.
const StoreType STORE_TYPE = precompute;


class CorrelatedStore {
protected:
  const size_t batch_size;
  const int server_num;
  const int serverfd;
  const bool lazy;  // Lazy (efficient but insecure) for behavior testing.

public:
  const bool do_fork;

  CorrelatedStore(
      const int serverfd, const int server_num, const size_t batch_size,
      const bool lazy = false, const bool do_fork = true)
  : batch_size(batch_size)
  , server_num(server_num)
  , serverfd(serverfd)
  , lazy(lazy)
  , do_fork(do_fork)
  {
    if (lazy) {
      std::cout << "Doing fast but insecure dabit computation" << std::endl;
    }
  }

  // Convert N single-bit values
  virtual fmpz_t* const b2a_single(const size_t N, const bool* const x) = 0;
  // Convert N values, with each x[i] having num_bits[i] bits.
  virtual fmpz_t* const b2a_multi(
      const size_t N, const size_t* const num_bits, const fmpz_t* const x) = 0;
  virtual ~CorrelatedStore() {};
};

class DaBitStore : public CorrelatedStore {
protected:
  std::queue<const DaBit* const> dabit_store;

public:
  DaBitStore(const int serverfd, const int server_num, const size_t batch_size,
      const bool lazy = false, const bool do_fork = true)
  : CorrelatedStore(serverfd, server_num, batch_size, lazy, do_fork)
  {
  }

  virtual ~DaBitStore() {
    while (!dabit_store.empty()) {
      const DaBit* const bit = dabit_store.front();
      dabit_store.pop();
      delete bit;
    }
  }

  // x, ret is [N]
  // Turns binary share x[i] into arith share ret[i]
  // Single bit. One round.
  fmpz_t* const b2a_single(const size_t N, const bool* const x);
  // Multiple bits. Use a dabit per bit in parallel, so one round
  fmpz_t* const b2a_multi(const size_t N, const size_t* const num_bits,
                          const fmpz_t* const x);

  // Check if enough to make n. Else, make or complain.
  // For batch generation to do it all at once.
  // Dependent on implementation / generation method.
  virtual void checkDaBits(const size_t n = 0) = 0;
  // Get a dabit from store. Checks for 1 first.
  const DaBit* const getDaBit();
};

// A Cache of correlated bits of different types
// Makes batch_size at once, when running low
class PrecomputeStore : public DaBitStore {
  // return N new daBits
  const DaBit* const * const generateDaBit(const size_t N);

  // add to the store.
  // Adds at least batch_size (or bool_batch_size), or n if bigger
  void addDaBits(const size_t n = 0);

  OT_Wrapper* const ot0;
  [[maybe_unused]] OT_Wrapper* const ot1;

public:

  PrecomputeStore(const int serverfd, const int server_num,
                  OT_Wrapper* const ot0, OT_Wrapper* const ot1,
                  const size_t batch_size,
                  const bool lazy = false, const bool do_fork = true)
  : DaBitStore(serverfd, server_num, batch_size, lazy, do_fork)
  , ot0(ot0)
  , ot1(ot1)
  {
  }

  ~PrecomputeStore() {};

  void printSizes();
  void maybeUpdate(); // Precompute if not enough.
  void checkDaBits(const size_t n = 0);
};

class OTCorrelatedStore : public CorrelatedStore {
  OT_Wrapper* const ot0;
  [[maybe_unused]] OT_Wrapper* const ot1;

public:

  OTCorrelatedStore(const int serverfd, const int server_num,
      OT_Wrapper* const ot0, OT_Wrapper* const ot1,
      const size_t batch_size = 0,
      const bool lazy = false, const bool do_fork = true)
  : CorrelatedStore(serverfd, server_num, batch_size, lazy, do_fork)
  , ot0(ot0)
  , ot1(ot1)
  {
  }

  ~OTCorrelatedStore() {};

  fmpz_t* const b2a_single(const size_t N, const bool* const x);
  fmpz_t* const b2a_multi(const size_t N, const size_t* const num_bits,
                          const fmpz_t* const x);

  fmpz_t* const b2a_ot(const size_t num_shares, const size_t num_values,
      const size_t* const num_bits,
      const fmpz_t* const x, const size_t mod);
};

// "legacy" EdaBit boolean triples.
// No longer relevant, as rounds scaled with num bits.
// Kept as legacy, since bool triple logic can be helpful
class BoolStore : public CorrelatedStore {
  std::queue<const BooleanBeaverTriple* const> btriple_store;
  void addBoolTriples(const size_t n = 0);

   OT_Wrapper* const ot0;
  OT_Wrapper* const ot1;

public:
  BoolStore(const int serverfd, const int server_num,
      OT_Wrapper* const ot0, OT_Wrapper* const ot1,
      const size_t batch_size,
      const bool lazy = false, const bool do_fork = true)
  : CorrelatedStore(serverfd, server_num, batch_size, lazy, do_fork)
  , ot0(ot0)
  , ot1(ot1)
  {
  }

  ~BoolStore();

  const BooleanBeaverTriple* const getBoolTriple();
  void checkBoolTriples(const size_t n = 0);

  void printSizes();
  void maybeUpdate();

  // x, y, ret are [N]
  // does ret[i] = x[i] * y[i], as shares
  const bool* const multiplyBoolShares(
      const size_t N, const bool* const x, const bool* const y);

  // x, y, z are [N][num_bits], ret is [N]
  // Treats x[i], y[i], z[i] as array of bits
  // sets z[i] and ret[i] as x[i] + y[i] and carry[i], as shares
  // Currently unused
  const bool* const addBinaryShares(
      const size_t N, const size_t* const num_bits,
      const bool* const * const x, const bool* const * const y,
      bool* const * const z);

  // Unimplemented
  fmpz_t* const b2a_single(const size_t N, const bool* const x);
  // Multiple bits. Use a dabit per bit in parallel, so one round
  fmpz_t* const b2a_multi(const size_t N, const size_t* const num_bits,
                          const fmpz_t* const x);
};


#endif
