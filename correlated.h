#ifndef CORRELATED_H
#define CORRELATED_H

/*
Correlated precomputes

CorrelatedStore object, which has synced correlated precomputes between servers
Also includes io objects for OT, which have their own precomputes

CorrelatedStore maintains a cache of precomputes, made in batch_size chunks at a time to reduce rounds
New batches are build either as it runs out, or by calling maybeUpdate

For now, only uses DaBits for b2a Share conversion.
boolean beaver triples are supported as they are straightforward, but not currently made.

*/

#include <emp-ot/emp-ot.h>
#include <emp-tool/emp-tool.h>
#include <queue>

#include "constants.h"
#include "ot.h"
#include "share.h"

// A Cache of correlated bits of different types
// Makes batch_size at once, when running low
class CorrelatedStore {
  const size_t batch_size;  // How many to make at once
  const int server_num;
  const int serverfd;

  // If lazy, does fast but insecure offline.
  const bool lazy;

  std::queue<const DaBit*> dabit_store;
  std::queue<const BooleanBeaverTriple*> btriple_store;

  // return N new daBits
  const DaBit* const * const generateDaBit(const size_t N);

  // add to the store.
  // Adds at least batch_size (or bool_batch_size), or n if bigger
  void addBoolTriples(const size_t n = 0);
  void addDaBits(const size_t n = 0);

  OT_Wrapper* const ot0;
  OT_Wrapper* const ot1;

public:

  CorrelatedStore(const int serverfd, const int idx,
                  OT_Wrapper* const ot0, OT_Wrapper* const ot1,
                  const size_t batch_size,
                  const bool lazy = false)
  : batch_size(batch_size)
  , server_num(idx)
  , serverfd(serverfd)
  , lazy(lazy)
  , ot0(ot0)
  , ot1(ot1)
  {
    if (lazy) {
      std::cout << "Doing fast but insecure dabit precomputes." << std::endl;
    }
  }

  ~CorrelatedStore();

  // get from store, and maybe add if necessary
  const BooleanBeaverTriple* const getBoolTriple();
  const DaBit* const getDaBit();

  void printSizes() const;
  // Precompute if not enough.
  void maybeUpdate();

  // check if enough to make n. if not, call add
  void checkBoolTriples(const size_t n = 0);
  void checkDaBits(const size_t n = 0);

  // compute with store elements. Does batches of size N.

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

  // x, ret is [N]
  // Turns binary share x[i] into arith share ret[i]
  // Single bit. One round.
  fmpz_t* const b2a_daBit_single(const size_t N, const bool* const x);
  // Multiple bits. Use a dabit per bit in parallel, so one round
  fmpz_t* const b2a_daBit_multi(const size_t N, const size_t* const num_bits,
                                const fmpz_t* const x);

  // Using intsum_ot, multiple bits
  // TODO: shift mod to fmpz
  fmpz_t* const b2a_ot(
      const size_t num_shares, const size_t num_values,
      const size_t* const num_bits, const fmpz_t* const shares,
      const size_t mod = 0);
};

#endif
