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
#include "proto.h"
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

public:

  // for OTs between servers
  NetIO* io0;
  NetIO* io1;

  // True to work correctly on large batches, and faster
  // False to do debugging
  const bool do_fork;

  CorrelatedStore(const int serverfd, const int idx,
                  const char* const server0_ip, const char* const server1_ip,
                  const size_t nbits, const size_t batch_size = 64,
                  const bool lazy = false, const bool do_fork = true)
  : batch_size(batch_size)
  , server_num(idx)
  , serverfd(serverfd)
  , nbits(nbits)
  , lazy(lazy)
  , bool_batch_size(batch_size * nbits)
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
    io0 = new NetIO(server_num == 0 ? nullptr : server0_ip, 60052, true);
    io1 = new NetIO(server_num == 1 ? nullptr : server1_ip, 60053, true);
  }

  ~CorrelatedStore();

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

  // get from store, and maybe add if necessary
  BooleanBeaverTriple* getBoolTriple();
  BeaverTriple* getTriple();
  DaBit* getDaBit();
  EdaBit* getEdaBit(const size_t num_bits);

  // Precompute if not enough
  void maybeUpdate();

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
  fmpz_t* b2a_daBit(const size_t N, const bool* const x);
  fmpz_t* b2a_edaBit(const size_t N, const size_t* const num_bits,
                     const fmpz_t* const x);

  // Unused
  // x2, xp, ret are [N]
  // ret[i] is if xor shares x2[i] and additive shares xp[i] are the same value
  // bool* validateSharesMatch(const size_t N, const size_t* const num_bits,
  //                           const fmpz_t* const x2, const fmpz_t* const xp);
};

#endif
