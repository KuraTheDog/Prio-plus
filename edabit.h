#ifndef EDABIT_H
#define EDABIT_H

#include <emp-ot/emp-ot.h>
#include <emp-tool/emp-tool.h>
#include <queue>

#include "proto.h"
#include "share.h"

// Return [z] = [x] and [y], consuming a boolean triple
bool multiplyBoolShares(const int serverfd, const int server_num, const bool x,
                        const bool y, BooleanBeaverTriple* const triple);

// Set [z] = [x] * [y], consuming an arithemtic triple
void multiplyArithmeticShares(const int serverfd, const int server_num, 
                              const fmpz_t x, const fmpz_t y, fmpz_t z, 
                              const BeaverTriple* const triple);

// Add two n-bit bool arrays x and y, fill in n bit array z
// Uses n boolean triples, and returns carry bit
bool addBinaryShares(const int serverfd, const int server_num, const size_t n,
                     const bool* const x, const bool* const y, bool* const z,
                     std::queue<BooleanBeaverTriple*> triples);

// Convert bit share [b]_2 into [b]_p using a daBit
void b2a_daBit(const int serverfd, const int server_num,
               const DaBit* const dabit, const bool x, fmpz_t& xp);

// Convert int share [x]_2 into [x]_p
// Assumes x, r < p/2 + 1
void b2a_edaBit(const int serverfd, const int server_num,
                const fmpz_t x, fmpz_t& xp, const EdaBit* const edabit,
                std::queue<BooleanBeaverTriple*> triples);

// Create a dabit, consuming an arithmetic triple
DaBit* generateDaBit(const int serverfd, const int server_num,
                     const BeaverTriple* const triple);

// Create a size n edaBit, consuming n boolean triples and a dabit
EdaBit* generateEdaBit(const int serverfd, const int server_num,const size_t n,
                       std::queue<BooleanBeaverTriple*> triples,
                       const DaBit* const dabit);

// Return if [x2]_2 and [xp]_p encode the same value, using an edaBit and n triples
bool validate_shares_match(const int serverfd, const int server_num,
                           const fmpz_t x2, const fmpz_t xp, const size_t n,
                           const EdaBit* const edabit,
                           std::queue<BooleanBeaverTriple*> triples);

// A Cache of correlated bits of different types
// Makes batch_size at once, when running low
class CorrelatedStore {
  const size_t batch_size;  // How many to make at once
  const int server_num;
  const int serverfd;
  const size_t num_bits;    // for size of edabits

  // If lazy, does fast but insecure offline.
  const bool lazy;

  // Since we use these a lot more, make much bigger batches at once.
  // num_bits * batch_size
  const size_t bool_batch_size;

  std::queue<EdaBit*> edabits;
  std::queue<DaBit*> dabits;

  // Maybe have these also pointers, for style consistency?
  std::queue<BooleanBeaverTriple*> bool_triples;
  std::queue<BeaverTriple*> triples;

  // for OTs between servers
  NetIO* io0;
  NetIO* io1;

public:

  CorrelatedStore(const int serverfd, const int idx, const char* const server0_ip, const char* const server1_ip, const size_t num_bits, const size_t batch_size = 64, const bool lazy = false) 
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

  ~CorrelatedStore() {
    while (!edabits.empty()) {
      EdaBit* bit = edabits.front();
      edabits.pop();
      delete bit;
    }
    while (!dabits.empty()) {
      DaBit* bit = dabits.front();
      dabits.pop();
      delete bit;
    }
    while (!triples.empty()) {
      BeaverTriple* triple = triples.front();
      triples.pop();
      delete triple;
    }
    while (!bool_triples.empty()) {
      BooleanBeaverTriple* triple = bool_triples.front();
      bool_triples.pop();
      delete triple;
    }
    delete io0;
    delete io1;
  }

  void addBoolTriples();
  void addTriples();
  void addDaBits();
  void addEdaBits();

  std::queue<BooleanBeaverTriple*> getBoolTriples(const size_t n);
  BeaverTriple* getTriple();
  DaBit* getDaBit();
  EdaBit* getEdaBit();

  // Precompute if not enough
  void maybeUpdate();
};

#endif
