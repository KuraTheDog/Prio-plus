#ifndef CORRELATED_H
#define CORRELATED_H

/* 
Correlated precomputes

CorrelatedStore object, which has synced correlated precomputes between servers
Also includes io objects for OT, which have their own precomputes

CorrelatedStore maintains a cache of precomputes, made in batch_size chunks at a time to reduce rounds
New batches are build either as it runs out, or by calling maybeUpdate

edaBit related logic for share conversion based on ia.cr/2020/338

Due to send buffers potentially filling up, it forks out a child to do sending, while parent receives
It also waits for the child to finish before exiting or moving to a substep that will send, to stay synced
*/

#include <emp-ot/emp-ot.h>
#include <emp-tool/emp-tool.h>
#include <queue>

#include "constants.h"
#include "he_triples.h"
#include "ot.h"
#include "share.h"

// A Cache of correlated bits of different types
// Makes batch_size at once, when running low
class CorrelatedStore {
  const size_t batch_size;  // How many to make at once
  const int server_num;
  const int serverfd;
  const size_t nbits;       // base # of bits. Makes this and *2

  // If lazy, does fast but insecure offline.
  const bool lazy;
  // If set, makes extra objects for precomputing. 
  // E.g. makes n edabits, then also n dabits in case the edabits run out
  // Allows for remaking edabits if they run out, without being from scratch
  const bool over_precompute;

  // Arithmetic triple generator
  ArithTripleGenerator* triple_gen = nullptr;

  // Since we use these a lot more, make much bigger batches at once.
  // nbits * batch_size
  const size_t bool_batch_size;

  std::queue<EdaBit*> edabit_store;    // nbits edabits
  std::queue<EdaBit*> edabit_store_2;  // 2 nbits edabits
  std::queue<DaBit*> dabit_store;
  std::queue<BooleanBeaverTriple*> btriple_store;
  std::queue<BeaverTriple*> atriple_store;

  // return N new daBits
  DaBit** generateDaBit(const size_t N);
  // return N new edaBits
  EdaBit** generateEdaBit(const size_t N, const size_t num_bits);

  // add to the store.
  // Adds at least batch_size (or bool_batch_size), or n if bigger
  void addBoolTriples(const size_t n = 0);
  void addTriples(const size_t n = 0);
  void addDaBits(const size_t n = 0);
  void addEdaBits(const size_t num_bits, const size_t n = 0);

  // check if enough to make n. if not, call add
  void checkBoolTriples(const size_t n = 0);
  void checkTriples(const size_t n = 0, const bool always = false);
  void checkDaBits(const size_t n = 0);
  void checkEdaBits(const size_t num_bits, const size_t n = 0);

  OT_Wrapper* const ot0;
  OT_Wrapper* const ot1;

public:

  // True to work correctly on large batches, and faster
  // False to do debugging
  const bool do_fork;

  CorrelatedStore(const int serverfd, const int idx,
                  OT_Wrapper* const ot0, OT_Wrapper* const ot1,
                  const size_t nbits, const size_t batch_size,
                  const bool lazy = false, const bool do_fork = true,
                  const bool over_precompute = true)
  : batch_size(batch_size)
  , server_num(idx)
  , serverfd(serverfd)
  , nbits(nbits)
  , lazy(lazy)
  , over_precompute(over_precompute)
  , bool_batch_size(batch_size * nbits)
  , ot0(ot0)
  , ot1(ot1)
  , do_fork(do_fork)
  {
    if (lazy) {
      std::cout << "Doing fast but insecure precomputes." << std::endl;
    } else if (fmpz_cmp_ui(Int_Modulus, 1ULL << 49) < 0) {
      std::cout << "Using PALISADE SHE arith triples" << std::endl;
      triple_gen = new ArithTripleGenerator(serverfd, server_num);
    } else {
      std::cout << "Mod big, using slower arith triples" << std::endl;
    }
  }

  ~CorrelatedStore();

  // get from store, and maybe add if necessary
  BooleanBeaverTriple* getBoolTriple();
  BeaverTriple* getTriple();
  DaBit* getDaBit();
  EdaBit* getEdaBit(const size_t num_bits);

  void printSizes();
  // Precompute if not enough.
  // if eda = true, generates edabits
  // otherwise makes dabits for dabit_multi
  void maybeUpdate(const bool using_eda = false);

  // compute with store elements. Does batches of size N.

  // x, y, ret are [N]
  // does ret[i] = x[i] * y[i], as shares
  bool* multiplyBoolShares(const size_t N,
                           const bool* const x, const bool* const y);
  fmpz_t* multiplyArithmeticShares(const size_t N,
                                   const fmpz_t* const x, const fmpz_t* const y);

  // x, y, z are [N][num_bits], ret is [N]
  // Treats x[i], y[i], z[i] as array of bits
  // sets z[i] and ret[i] as x[i] + y[i] and carry[i], as shares
  bool* addBinaryShares(const size_t N, const size_t* const num_bits,
                        const bool* const * const x, const bool* const * const y,
                        bool* const * const z);

  // x, ret is [N]
  // Turns binary share x[i] into arith share ret[i]
  // Single bit. One round.
  fmpz_t* b2a_daBit_single(const size_t N, const bool* const x);
  // Multiple bits. Use a dabit per bit in parallel, so one round
  fmpz_t* b2a_daBit_multi(const size_t N, const size_t* const num_bits,
                          const fmpz_t* const x);
  // Multiple bits. Uses one edabit, which requires L rounds for L bits.
  fmpz_t* b2a_edaBit(const size_t N, const size_t* const num_bits,
                     const fmpz_t* const x);

  // Using intsum_ot, multiple bits
  // TODO: shift mod to fmpz
  fmpz_t* b2a_ot(const size_t num_shares, const size_t num_values, 
                 const size_t* const num_bits, const fmpz_t* const shares,
                 const size_t mod = 0);
};

#endif
