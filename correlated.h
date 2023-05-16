#ifndef CORRELATED_H
#define CORRELATED_H

/*
Correlated precomputes
synced correlated precomputed values between servers
Maintaines a cache of precomputed correlated randomness.

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
Pure live OT based B2A, via intsum_ot
(Implicitly has OT extension based precompute)

ShareConverter
* General B2A share conversion skeleton
- CorrelatedStore
  * Dabit storage
  * triple storage (arith and bool)
  * Other triple-based computation
  - PrecomputeStore
    * generation via Offline precompute
    - validated
      * generation via client provided and server validated
      * still uses "offline" in case of bad values
- OT
  * OT based B2A only
*/

#include <emp-ot/emp-ot.h>
#include <emp-tool/emp-tool.h>
#include <queue>

#include "constants.h"
// #include "he_triples.h"
#include "ot.h"
#include "share.h"

// Selector
enum StoreType {
  precompute_store,
  ot_store,
  validate_store,
};

// Global for now. Can have client/server sync it or be a param file, but not necessary.
// const StoreType STORE_TYPE = precompute_store;
// const StoreType STORE_TYPE = ot_store;
const StoreType STORE_TYPE = validate_store;


class ShareConverter {
protected:
  const int server_num;
  const int serverfd;

public:

  ShareConverter(
      const int serverfd, const int server_num)
  : server_num(server_num)
  , serverfd(serverfd)
  {};

  // Convert N single-bit values
  virtual int b2a_single(const size_t N, const bool* const x, fmpz_t* const xp) = 0;
  // Convert N values, with each x[i] having num_bits[i] bits.
  virtual int b2a_multi(const size_t N, const size_t* const num_bits,
                        const fmpz_t* const x, fmpz_t* const xp) = 0;
  virtual ~ShareConverter() {};
};

class CorrelatedStore : public ShareConverter {
protected:
  OT_Wrapper* const ot0;
  [[maybe_unused]] OT_Wrapper* const ot1;

  std::queue<const DaBit*> dabit_store;
  std::queue<const BooleanBeaverTriple*> btriple_store;
  std::queue<const BeaverTriple*> atriple_store;

  // Get: Return 1 (check first, which may generate)
  virtual const DaBit* const get_DaBit();
  virtual const BeaverTriple* const get_Triple();
  virtual const BooleanBeaverTriple* const get_BoolTriple();

  // Currently only used for heavy eval, which is independent of # clients
  // Approx 41 * b * nbits_mod, where b is #bits of run
  // Here estimating a rounds of 16 bits (or 4 of 4, etc)
  // TODO: This will get deprecated with garbled
  const size_t triples_batch_size = 41 * 16 * nbits_mod;

public:
  CorrelatedStore(const int serverfd, const int server_num,
                  OT_Wrapper* const ot0, OT_Wrapper* const ot1)
  : ShareConverter(serverfd, server_num)
  , ot0(ot0)
  , ot1(ot1)
  {
  }

  virtual ~CorrelatedStore() {
    while (!dabit_store.empty()) {
      const DaBit* const bit = dabit_store.front();
      dabit_store.pop();
      delete bit;
    }
    while (!btriple_store.empty()) {
      const BooleanBeaverTriple* const triple = btriple_store.front();
      btriple_store.pop();
      delete triple;
    }
    while (!atriple_store.empty()) {
      const BeaverTriple* const triple = atriple_store.front();
      atriple_store.pop();
      delete triple;
    }
  }

  // Single bit. One round.
  int b2a_single(const size_t N, const bool* const x, fmpz_t* const xp);
  // Multiple bits. Use a dabit per bit in parallel, so one round
  int b2a_multi(const size_t N, const size_t* const num_bits,
                const fmpz_t* const x, fmpz_t* const xp);

  // x, y, z are [N]
  // does z[i] = x[i] * y[i], as shares
  int multiply_BoolShares(const size_t N,
      const bool* const x, const bool* const y, bool* const z);
  int multiply_ArithmeticShares(
      const size_t N, const fmpz_t* const x, const fmpz_t* const y,
      fmpz_t* const z);
  // N inputs, each B bits large
  // b, x, and z are [N * B]
  // Distinction of N vs B only comes in for the optional "valid" array.
  // If z_inv (size N*B): Does a second set multiplied by 1-b instead.
  //   maintains same rounds, just "larger" OT's
  // If valid, it's a [N] array for which inputs are valid
  int multiply_BoolArith(
      const size_t N, const size_t B, const bool* const b, const fmpz_t* const x,
      fmpz_t* const z, fmpz_t* const z_inv = nullptr, const bool* const valid = nullptr);
  // b is [B] rather than [N * B]
  // Grouped by B first, so [x0, ..., x(B-1)] for N=0
  // So [x0, xB, x2B, ...] multiplied by b0.
  // For B-wide masking, mainly
  int multiply_BoolArithFlat(
      const size_t N, const size_t B, const bool* const b_flat, const fmpz_t* const x,
      fmpz_t* const z, fmpz_t* const z_inv = nullptr, const bool* const valid = nullptr);

  // x, y, z are [N][num_bits], carry is [N]
  // Treats x[i], y[i], z[i] as array of bits
  // sets z[i] and carry[i] as x[i] + y[i] and carry[i], as shares
  // Currently unused
  [[maybe_unused]] int add_BinaryShares(
        const size_t N, const size_t* const num_bits,
        const bool* const * const x, const bool* const * const y,
        bool* const * const z, bool* const carry);

  // N inputs of b values
  // [xy]: x = which bucket is nonzero, nonzero is -1 if y = 1
  // 00 = (0, 1), 01 = (0, -1), 10 = (1, 0), 11 = (-1, 0)
  // |x| = |y| = N * b
  // Index with i * b + j: (x0, ..., xb-1) for number 0.
  // Makes copying in new shares easier.
  // N valid, if invalid just contributes 0
  int heavy_convert(const size_t N, const size_t b,
                    const bool* const x, const bool* const y,
                    const bool* const valid,
                    fmpz_t* const bucket0, fmpz_t* const bucket1);

  // Q "parallel" runs of heavy_convert
  // N inputs, Q copies, D depth, B substreams
  // |x| = |y| = N * Q * D
  // |mask| = N * Q * M: Substream select
  // Valid size N
  // buckets size : Q * M * D
  // Order: N, Q, M, D
  //   x,y = [over N (over Q (over D))]
  //   valid = over N
  //   mask  = [over N (over Q (over M))]
  //   buckets = [over Q (over M (over D))]
  int heavy_convert_mask(
      const size_t N, const size_t Q, const size_t M, const size_t D,
      const bool* const x, const fmpz_t* const y_p, const bool* const mask,
      const bool* const valid, fmpz_t* const bucket0, fmpz_t* const bucket1);

  /* Comparisons.
  // Assumes modulus is large enough, so <N/2 is positive, >N/2 is negative (x-N)
  Return true if first < second. (so is index of larger).
  Equality is merged in, since rare/should not happen case so far.
    Specifically, if x < y then 1, if x >= y then 0
  Returns additive shares of the comparsion
  */

  // Clear: Reveals info, not secure
  // True if [x] > c
  bool* cmp_c_clear(const size_t N, const fmpz_t* const x, const fmpz_t* const c);
  // True if [x] > [y]
  int cmp(const size_t N, const fmpz_t* const x, const fmpz_t* const y, fmpz_t* const ans);
  // Given [x], sets out = [|x|] via negating if [x] < 0
  int abs(const size_t N, const fmpz_t* const x, fmpz_t* const abs_x);
  // True if [|x|] > [|y|]
  int abs_cmp(const size_t N,
              const fmpz_t* const x, const fmpz_t* const y, fmpz_t* const ans);
  // x, y are Nxb shares of N total b-bit numbers.
  // I.e. x[i,j] is additive share of bit j of number i.
  // Returns additive shares of [x < y]
  // Uses 3 triples per bit to compare
  int cmp_bit(const size_t N, const size_t b,
                  const fmpz_t* const x, const fmpz_t* const y, fmpz_t* const ans);
  // Generate N random b-bit numbers r, with bitshares rB
  // b is the # of bits in the int modulus
  // r is N values [r], returns rB is N*b shares of bits
  // For each gen, uses a dabit to gen and a cmp_bit (3 triple) to check
  int gen_rand_bitshare(const size_t N, fmpz_t* const r, fmpz_t* const rB);
  // Extracts shares of the least significant bit of additive secret [x]
  int LSB(const size_t N, const fmpz_t* const x, fmpz_t* const x0);
  // Returns share of
  //    1 if x is negative (x > p/2)
  //    0 if x is positive (x < p/2)
  // x is N b-bit numbers as additive shares
  int is_negative(const size_t N, const fmpz_t* const x, fmpz_t* const ans);

  // Check: Check if enough to make n.
  // If not enough, generates enough (can over-generate)
  virtual void check_DaBits(const size_t n = 0) = 0;
  virtual void check_Triples(const size_t n = 0) = 0;
  virtual void check_BoolTriples(const size_t n = 0) = 0;
};

// A Cache of correlated bits of different types
// Makes batch_size at once, when running low
class PrecomputeStore : public CorrelatedStore {
protected:
  const size_t batch_size;
  const bool lazy;  // Lazy (efficient but insecure) generation for behavior testing.

  // Securely create N new correlated items
  const DaBit* const * const gen_DaBits(const size_t N);
  const DaBit* const * const gen_DaBits_lazy(const size_t N);

  // Add (securely or lazy generated) N items to the store, at least batch_size many
  // TODO: return sent bytes
  void add_DaBits(const size_t n = 0);
  void add_Triples(const size_t n = 0);
  void add_BoolTriples(const size_t n = 0);

public:

  PrecomputeStore(const int serverfd, const int server_num,
                  OT_Wrapper* const ot0, OT_Wrapper* const ot1,
                  const size_t batch_size,
                  const bool lazy = false)
  : CorrelatedStore(serverfd, server_num, ot0, ot1)
  , batch_size(batch_size)
  , lazy(lazy)
  {
    if (lazy) {
      std::cout << "Doing fast but insecure correlated value computation" << std::endl;
    }
  };

  ~PrecomputeStore() {};

  virtual void print_Sizes() const;
  virtual void maybe_Update(); // Precompute if not enough.

  virtual void check_DaBits(const size_t n = 0);
  virtual void check_Triples(const size_t n = 0);
  virtual void check_BoolTriples(const size_t n = 0);
};

class OTCorrelatedStore : public ShareConverter {
  OT_Wrapper* const ot0;
  [[maybe_unused]] OT_Wrapper* const ot1;

public:

  OTCorrelatedStore(const int serverfd, const int server_num,
      OT_Wrapper* const ot0, OT_Wrapper* const ot1)
  : ShareConverter(serverfd, server_num)
  , ot0(ot0)
  , ot1(ot1)
  {
  }

  ~OTCorrelatedStore() {};

  int b2a_single(const size_t N, const bool* const x, fmpz_t* const xp);
  // Multiple bits. Use a dabit per bit in parallel, so one round
  int b2a_multi(const size_t N, const size_t* const num_bits,
                const fmpz_t* const x, fmpz_t* const xp);

  int b2a_ot(const size_t num_shares, const size_t num_values,
      const size_t* const num_bits,
      const fmpz_t* const x, fmpz_t* const xp,
      const size_t mod = 0);
};

#endif
