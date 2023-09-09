#ifndef CORRELATED_H
#define CORRELATED_H

/*
Correlated precomputes
synced correlated precomputed values between servers
Maintaines a cache of precomputed correlated randomness.

PrecomputeStore:
Made in at least batch_size chunks at a time.
New batches are build either as it runs out, or by calling maybe_Update
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

/*
Stage atomicizing. So can do multiple send/rec in parallel
See inside main functions for uses. (e.g. inside b2a_single for how to use)
1) Setup params
2) Check
3) Declare variables. Things passed to send/rec. Need to be new/nullptr outside.
4) setup
5) send + rec
6+) process
7+) send + rec
...
n-1) finish(params, vars)
n) delete vars. could fit into finish, but clearer outside
*/

#include <emp-ot/emp-ot.h>
#include <emp-tool/emp-tool.h>
#include <queue>

#include "constants.h"
// #include "he_triples.h"
#include "net_share.h"
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

// Useful in destructors
// Clear a queue of pointers by deleting all of them.
template <class T>
void clear_queue(std::queue<T*> q){
  while(!q.empty()) {
    const T* const x = q.front();
    q.pop();
    delete x;
  }
};

// Helper: x is N copies len a, y is N copies len b.
// Fill out x_ext, and y_ext as N x a x b
// So that x_ext * y_ext = (x cross y)
// Could be template, but just bool used (since doesn't work on fmpz)
void cross_fill_bool(
    const size_t N, const size_t a, const size_t b,
    const bool* const x, const bool* const y,
    bool* const x_ext, bool* const y_ext);


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
  virtual int64_t b2a_single(const size_t N, const bool* const x, fmpz_t* const xp) = 0;
  // Convert N values, with each x[i] having num_bits bits.
  virtual int64_t b2a_multi(const size_t N, const size_t num_bits,
                            const uint64_t* const x, fmpz_t* const xp) = 0;
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

  // Helper, not expected to be used outside.
  /* N inputs across m accumulators
    If mask: values is (m, xm, ym, xym), len 4NM
    If no mask, values is (xm, ym, xym), len 3NM
    Valid length N, Buckets each length M
  */
public:
  CorrelatedStore(const int serverfd, const int server_num,
                  OT_Wrapper* const ot0, OT_Wrapper* const ot1)
  : ShareConverter(serverfd, server_num)
  , ot0(ot0)
  , ot1(ot1)
  {
  }

  virtual ~CorrelatedStore() {
    clear_queue(dabit_store);
    clear_queue(btriple_store);
    clear_queue(atriple_store);
  }

  // Single bit value.
  // One round, for each consuming a dabit and sending 1 bit.
  void b2a_single_setup(
      const size_t N, const bool* const x, fmpz_t* const xp, bool* const v);
  void b2a_single_finish(
      const size_t N, fmpz_t* const xp, const bool* const v, const bool* const v_other);
  int64_t b2a_single(const size_t N, const bool* const x, fmpz_t* const xp);

  // Multiple bits values.
  // One round, for each bit consuming a dabit and sending 1 bit.
  // x[i] is length num_bits boolean share
  // Different num_bits can be done with parallel runs.
  void b2a_multi_setup(
      const size_t N, const size_t num_bits,
      const uint64_t* const x, fmpz_t* const flat_xp, bool* const v);
  void b2a_multi_finish(
      const size_t N, const size_t num_bits,
      fmpz_t* const xp, fmpz_t* const flat_xp, const bool* const v, const bool* const v_other);
  int64_t b2a_multi(const size_t N, const size_t num_bits,
      const uint64_t* const x, fmpz_t* const xp);

  // x, y, z are [N]
  // does z[i] = x[i] * y[i], as shares
  // Note: de size 2 * N
  /* TODO: It would be cool to have a "binary-fold" multiply list
     E.g. for 5 to 8 items, do in 3 rounds via pairing
     For "heavy mask" though, would require weirdness with unfolding
     since "variable" number of rounds.
  */
  void multiply_BoolShares_setup(const size_t N,
      const bool* const x, const bool* const y, bool* const z, bool* const de);
  void multiply_BoolShares_finish(const size_t N,
      const bool* const x, const bool* const y, bool* const z,
      const bool* const de, const bool* const de_other);
  int64_t multiply_BoolShares(const size_t N,
      const bool* const x, const bool* const y, bool* const z);
  int64_t multiply_ArithmeticShares(const size_t N,
      const fmpz_t* const x, const fmpz_t* const y, fmpz_t* const z);
  // x is [N x a], y is [M x b], cross multiply into z as [N * (a * b)]
  // Each of N inputs gives (a x b) cross
  // x gives column, y is row
  // so z = [x[0] * y, x[1] * y, ...] for input 1, then same for input 2, etc.
  int64_t multiply_BoolShares_cross(
      const size_t N, const size_t a, const size_t b,
      const bool* x, const bool* y, bool* const z);
  // N inputs, each B bits large
  // b, x, and z are [N * B]
  // Distinction of N vs B only comes in for the optional "valid" array.
  // If optional z_inv (size N*B): Does a second set multiplied by 1-b instead.
  //   maintains same rounds/number of OTs, just "larger" OT's (extra params)
  // If valid, it's a [N] array for which inputs are valid
  // Uses N * B OTs, no correlated
  int64_t multiply_BoolArith(
      const size_t N, const size_t B, const bool* const b, const fmpz_t* const x,
      fmpz_t* const z, fmpz_t* const z_inv = nullptr, const bool* const valid = nullptr);
  // b is [B] rather than [N * B]
  // Grouped by B first, so [x0, ..., x(B-1)] for N=0
  // So [x0, xB, x2B, ...] multiplied by b0.
  // For B-wide masking, mainly
  int64_t multiply_BoolArithFlat(
      const size_t N, const size_t B, const bool* const b_flat, const fmpz_t* const x,
      fmpz_t* const z, fmpz_t* const z_inv = nullptr, const bool* const valid = nullptr);

  // x, y, z are [N][num_bits], carry is [N]
  // Treats x[i], y[i], z[i] as array of bits
  // sets z[i] and carry[i] as x[i] + y[i] and carry[i], as shares
  // Currently unused
  [[maybe_unused]] int64_t add_BinaryShares(const size_t N, const size_t* const num_bits,
      const bool* const * const x, const bool* const * const y,
      bool* const * const z, bool* const carry);

  void heavy_accumulate(const size_t N, const size_t M,
      const fmpz_t* const values, const bool* const valid,
      fmpz_t* const bucket0, fmpz_t* const bucket1, const bool use_mask = true);

  /* N inputs (each may be valid or not) of b values
    [xy]: x = which bucket is nonzero, nonzero is -1 if y = 1
    00 = (0, 1), 01 = (0, -1), 10 = (1, 0), 11 = (-1, 0)
    |x| = |y| = N * b
    Index with i * b + j: (x0, ..., xb-1) for number 0.
    Accumulates into the buckets sized b.
  */
  // One round bool mult (N*b), one round of B2A (3*N*b)
  int64_t heavy_convert(const size_t N, const size_t b,
      const bool* const x, const bool* const y, const bool* const valid,
      fmpz_t* const bucket0, fmpz_t* const bucket1);
  // One round of B2A (N*b), one round of OTs (N*b)
  int64_t heavy_convert_ot(const size_t N, const size_t b,
      const bool* const x, const bool* const y, const bool* const valid,
      fmpz_t* const bucket0, fmpz_t* const bucket1);

  /* Mask versions
    For each of N inputs,
      values in groups sized D
        [xy]: x = which bucket is nonzero, nonzero is -1 if y = 1
      cross multiply against a mask sized M
      repeat Q times
      accumulate across Q * M * D
    Note that it repeats QN times, distinction is only for accumulation

    Sizes (and ordering)
    x, y = N (Q (D))
    valid = N
    mask = N (Q (M))
    buckets = Q (M (D))
    intermediate N (Q (M (D)))

    Note: Could fold in "valid" short circuits more.
    Won't affect runtime much since done at end anyways and branch prediction.
    Less memory copying if invalid.
  */
  // OT version.
  // Requires one round of B2A (NQD) outside: (1-2y)
  // Two rounds of Ots (NQMD each): first * mask, second a double * (x, 1-x)
  int64_t heavy_convert_mask(
      const size_t N, const size_t Q, const size_t M, const size_t D,
      const bool* const x, const fmpz_t* const y_p, const bool* const mask,
      const bool* const valid, fmpz_t* const bucket0, fmpz_t* const bucket1);
  /*
    Specific ones for levels of binary tree
    Due to binary tree packing of bool mults (e.g. 8 done in 3 rather than 7)
    for some "in between" amounts, can do mults outside
  */
  /* Single mask
    Two round bool mult (3NQMD), one round b2a (4 NQMD)
    First mult to get xm, ym. Second mult to get xym
    Then B2A on m, xm, ym, xym
    Math: (1-2y)xm = xm - 2xym, and (1-2y)(1-x)m = m-xm-2ym-2xym
    VS OT: More b2a (extra 4M factor).
      but then cheaper bool mult (factor 3), versus expensive OT (factor 2)
  */
  int64_t heavy_convert_mask_one(
      const size_t N, const size_t Q, const size_t M, const size_t D,
      const bool* const x, const bool* const y, const bool* const mask,
      const bool* const valid, fmpz_t* const bucket0, fmpz_t* const bucket1);
  /* Double mask
    If mult masks first then feed into above, is +1 mult round, for 3 mults
    Instead can binary fold rounds to have less.
    This is the core payoff.
    (1-2y)xm1m2 = xm1m2 - 2xym1m2
    (1-2y)(1-x)m1m2 = m1m2 - xm1m2 - 2ym1m2 + 2xym1m2
    Mult 1: m1m2, xy (bonus: not full dimensions, extend after)
    Mult 2: x(lm), y(lm), (xy)(lm)
    One less round:
      OT: Cross (N Q M1 M2), b2A (N Q D), OT (N Q M1 M2 D), OT (N Q M1 M2 D)
      This: Mult (N Q (D + M1 M2)), Mult (3 N Q M1 M2 D), b2A (4 N Q M1 M2 D)

    Note:
      Here, mask1 has size N (Q (M1)), and mask2 is N (Q (M2))
      Together, the combined mask is N Q M1 M2
      In use case however, mask 2 is just N M2, no Q.
      But since everything else has Q somewhere, it'll always be added in cross multiply.
      Hence, there's no benefit to the code supporting the smaller mask2.

    M2 fitter for N M2 -> N Q M2
    for (unsigned int i = 0; i < N * Q; i++)
      memcpy(&mask2_ext[i * M2], &mask2[(i / Q) * M2], M2);
  */
  int64_t heavy_convert_mask_two(
      const size_t N, const size_t Q, const size_t M1, const size_t M2, const size_t D,
      const bool* const x, const bool* const y, const bool* const mask1, const bool* const mask2,
      const bool* const valid, fmpz_t* const bucket0, fmpz_t* const bucket1);

  /* Comparisons.
  Assumes modulus is large enough, so <N/2 is positive, >N/2 is negative (x-N)
  Return true if first < second. (so is index of larger).
  Equality is merged in, since rare/should not happen case so far.
    Specifically, if x < y then 1, if x >= y then 0
  Returns additive shares of the comparsion.

  Currently mainly replaced by garble (eval_heavy).
  */

  // Clear: Reveals info, not secure
  // True if [x] > c
  bool* cmp_c_clear(const size_t N, const fmpz_t* const x, const fmpz_t* const c);
  // True if [x] > [y]
  int64_t cmp(const size_t N, const fmpz_t* const x, const fmpz_t* const y, fmpz_t* const ans);
  // Given [x], sets out = [|x|] via negating if [x] < 0
  int64_t abs(const size_t N, const fmpz_t* const x, fmpz_t* const abs_x);
  // True if [|x|] > [|y|]
  int64_t abs_cmp(const size_t N,
              const fmpz_t* const x, const fmpz_t* const y, fmpz_t* const ans);
  // x, y are Nxb shares of the bits of N total b-bit numbers.
  // I.e. x[i,j] is additive share of bit j of number i.
  // Returns additive shares of [x < y] (NOTE: opposite of cmp, for convenience)
  // Uses 3 triples per bit to compare
  int64_t cmp_bit(const size_t N, const size_t b,
              const fmpz_t* const x, const fmpz_t* const y, fmpz_t* const ans);
  /* Generate N random b-bit numbers r, with bits rB
    Both are additive shares
    b is the # of bits in the int modulus
    For each n, makes 1 share r, and b shares rB_i of the bits of r
    r = sum_i rB_i * 2^i
    For each gen, uses a dabit to gen and a cmp_bit (3 triple) to check
    TODO: these can be pre-computed.
    Similar to edabit, but not quite.
    Bit shares here are additive, while edabit bit shares are boolean.
  */
  int64_t gen_rand_bitshare(const size_t N, fmpz_t* const r, fmpz_t* const rB);
  // Extracts shares of the least significant bit of additive secret [x]
  int64_t LSB(const size_t N, const fmpz_t* const x, fmpz_t* const x0);
  // Returns share of
  //    1 if x is negative (x > p/2)
  //    0 if x is positive (x < p/2)
  // x is N b-bit numbers as additive shares
  int64_t is_negative(const size_t N, const fmpz_t* const x, fmpz_t* const ans);

  // Check: Check if enough to make n.
  // If not enough, generates enough (can over-generate)
  virtual int64_t check_DaBits(const size_t n = 0) = 0;
  virtual void check_Triples(const size_t n = 0) = 0;
  virtual void check_BoolTriples(const size_t n = 0) = 0;

  virtual void print_Sizes() const = 0;
};

// A Cache of correlated bits of different types
// Makes batch_size at once, when running low
class PrecomputeStore : public CorrelatedStore {
protected:
  const size_t batch_size;
  const bool lazy;  // Lazy (efficient but insecure) generation for behavior testing.
  mutable flint_rand_t lazy_seed;  // Since seeds are modified when used

  // Securely create N new correlated items
  // Pass in pointer to array. Methods make new.
  int64_t gen_DaBits(const size_t N, DaBit** const dabit);
  int64_t gen_DaBits_lazy(const size_t N, DaBit** const dabit) const;

  int64_t gen_BoolTriple_lazy(const size_t N, BooleanBeaverTriple** const triples) const;

  // Add (securely or lazy generated) N items to the store, at least batch_size many
  // TODO: return sent bytes
  int64_t add_DaBits(const size_t n = 0);
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
      std::cout << "Doing fast but insecure correlated value generation" << std::endl;
      flint_randinit(lazy_seed);
      if (server_num == 0) {
        send_seed(serverfd, lazy_seed);
      } else {
        recv_seed(serverfd, lazy_seed);
      }
    }
  };

  ~PrecomputeStore() {
    if (lazy) flint_randclear(lazy_seed);
  };

  virtual void print_Sizes() const;
  virtual void maybe_Update(); // Precompute if not enough.

  virtual int64_t check_DaBits(const size_t n = 0);
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

  int64_t b2a_single(const size_t N, const bool* const x, fmpz_t* const xp);
  // Multiple bits. Use a dabit per bit in parallel, so one round
  int64_t b2a_multi(const size_t N, const size_t num_bits,
      const uint64_t* const x, fmpz_t* const xp);

  // Num_bits left as array, since harder to unpack OT
  int64_t b2a_ot(const size_t num_shares, const size_t num_values,
      const size_t* const num_bits,
      const uint64_t* const x, fmpz_t* const xp,
      const size_t mod = 0);

  void print_Sizes() const {};
  void maybe_Update() {};
};

#endif
