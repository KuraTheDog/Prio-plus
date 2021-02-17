#ifndef EDABIT_H
#define EDABIT_H

#include <emp-ot/emp-ot.h>
#include <emp-tool/emp-tool.h>
#include <queue>

#include "proto.h"
#include "share.h"

// A Cache of correlated bits of different types
// Makes batch_size at once, when running low
class CorrelatedStore {
// public:
  const size_t batch_size;  // How many to make at once
  const int server_num;
  const int serverfd;
  const size_t num_bits;    // for size of edabits

  // If lazy, does fast but insecure offline.
  const bool lazy;

  // Since we use these a lot more, make much bigger batches at once.
  // num_bits * batch_size
  const size_t bool_batch_size;

  std::queue<EdaBit*> edabit_store;
  std::queue<DaBit*> dabit_store;
  std::queue<BooleanBeaverTriple*> btriple_store;
  std::queue<BeaverTriple*> atriple_store;

  // for OTs between servers
  NetIO* io0;
  NetIO* io1;

public:

  CorrelatedStore(const int serverfd, const int idx,
                  const char* const server0_ip, const char* const server1_ip,
                  const size_t num_bits, const size_t batch_size = 64,
                  const bool lazy = false) 
  : batch_size(batch_size)
  , server_num(idx)
  , serverfd(serverfd)
  , num_bits(num_bits)
  , lazy(lazy)
  , bool_batch_size(2 * batch_size * num_bits)
  {
    if (lazy)
      std::cout << "Doing fast but insecure precomputes." << std::endl;
    io0 = new NetIO(server_num == 0 ? nullptr : server0_ip, 60051, true);
    io1 = new NetIO(server_num == 1 ? nullptr : server1_ip, 60052, true);
  }

  ~CorrelatedStore();

  // return N new daBits
  DaBit** generateDaBit(const size_t N);
  // return N new edaBits
  EdaBit** generateEdaBit(const size_t N);

  // add to the store.
  // Adds at least batch_size (or bool_batch_size), or n if bigger
  void addBoolTriples(const size_t n = 0);
  void addTriples(const size_t n = 0);
  void addDaBits(const size_t n = 0);
  void addEdaBits(const size_t n = 0);

  // get from store, and maybe add if necessary
  BooleanBeaverTriple* getBoolTriple();
  BeaverTriple* getTriple();
  DaBit* getDaBit();
  EdaBit* getEdaBit();

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
  bool* addBinaryShares(const size_t N,
                        const bool* const * const x, const bool* const * const y,
                        bool* const * const z);

  // x, ret is [N]
  // Turns binary share x[i] into arith share ret[i]
  fmpz_t* b2a_daBit(const size_t N, const bool* const x);
  fmpz_t* b2a_edaBit(const size_t N, const fmpz_t* const x);

  // x2, xp, ret are [N]
  // ret[i] is if xor shares x2[i] and additive shares xp[i] are the same value
  bool* validateSharesMatch(const size_t N,
                            const fmpz_t* const x2, const fmpz_t* const xp);
};

#endif
