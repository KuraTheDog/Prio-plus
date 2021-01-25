#ifndef EDABIT_H
#define EDABIT_H

#include <queue>

#include "share.h"

// Return [z] = [x] and [y], consuming a boolean triple
bool multiplyBoolShares(const int serverfd, const int server_num, const bool x, const bool y, const BooleanBeaverTriple triple);

// Set [z] = [x] * [y], consuming an arithemtic triple
void multiplyArithmeticShares(const int serverfd, const int server_num, const fmpz_t x, const fmpz_t y, fmpz_t z, const BeaverTriple* const triple);

// Add two n-bit bool arrays x and y, fill in n bit array z
// Uses n triples, and returns carry bit
bool addBinaryShares(const int serverfd, const int server_num, const size_t n, const bool* const x, const bool* const y, bool* const z, const BooleanBeaverTriple* const triples);

// Convert bit share [b]_2 into [b]_p using a daBit
void b2a_daBit(const int serverfd, const int server_num, const DaBit* const dabit, const bool x, fmpz_t& xp);

// Convert int share [x]_2 into [x]_p
// Assumes x, r < p/2 + 1
void b2a_edaBit(const int serverfd, const int server_num, const EdaBit* const edabit, const fmpz_t x, fmpz_t& xp, const BooleanBeaverTriple* const triples);

// Create a dabit, consuming an arithmetic triple
DaBit* generateDaBit(const int serverfd, const int server_num, const BeaverTriple* const triple);

// Create a size n edaBit, consuming n boolean triples and a dabit
EdaBit* generateEdaBit(const int serverfd, const int server_num, const size_t n, const BooleanBeaverTriple* const triples, const DaBit* const dabit);

// Return if [x2]_2 and [xp]_p encode the same value, using an edaBit and n triples
bool validate_shares_match(const int serverfd, const int server_num, const fmpz_t x2, const fmpz_t xp, const size_t n, const EdaBit* const edabit, const BooleanBeaverTriple* const triples);

// A Cache of correlated bits of different types
// Makes batch_size at once, when running low
struct CorrelatedStore {
  const size_t batch_size;  // How many to make at once
  const int server_num;
  const int serverfd;
  const size_t num_bits;    // for size of edabits

  std::queue<EdaBit*> edabits;
  std::queue<DaBit*> dabits;

  // Maybe have these also pointers, for style consistency?
  std::queue<BooleanBeaverTriple> bool_triples;
  std::queue<BeaverTriple*> triples;

  CorrelatedStore(const int serverfd, const int idx, const size_t num_bits, const size_t batch_size = 64) 
  : batch_size(batch_size)
  , server_num(idx)
  , serverfd(serverfd)
  , num_bits(num_bits)
  {}

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
  }

  void addBoolTriples();
  void addTriples();
  void addDaBits();
  void addEdaBits();

  BooleanBeaverTriple* getBoolTriples(const size_t n);
  BeaverTriple* getTriple();
  DaBit* getDaBit();
  EdaBit* getEdaBit();
};

#endif
