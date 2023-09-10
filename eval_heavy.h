#ifndef EVAL_HEAVY_H
#define EVAL_HEAVY_H

#include "constants.h"
#include "hash.h"
#include "heavy.h"
#include "utils.h"

#include "emp-sh2pc/emp-sh2pc.h"

/*
Uses MPC circuit to evaluate heavy
Evaluates on Count-min

For the garble wrapping:
One time start:
- Set up NetIO (e.g. ot0->io)
- setup_semi_honest(garbleIO, server_num + 1);
Run:
- full_heavy_extract(...)
- garbleIO->flush()
End:
finalize_semi_honest();
*/

/*
Pieces can be abstracted out
But for now, in one struct so easier to access various attributes.
Mostly one-time, but can call with multiple K I guess.
Can theoretically be one big function with sub-functions.
Currently split for clarity of use / debugging.

There may be better ways of de-duping, but this should be fine.
- Challenge is to keep it oblivious. same work if all dupes and no dupes.
- Can't track how many distinct values seen
- Can't pre-prune list of candidates.

Note: emp::ALICE is 1, emp::BOB is 2. (Using server_num + 1)
*/

/* So that we can sort by (freq, value)
  Comparison returns Bit, since that's what Integer does
  And just needs geq, since that builds operators in template
  Operators then to support sorting

  So that same values are next to each other if there's shared frequency when sorting

  Inherit  : public Swappable<IntegerPair>, Comparable<IntegerPair> ?
*/
class IntegerPair{
public:
  Integer x;
  Integer y;

  IntegerPair(){};

  IntegerPair(Integer x, Integer y)
  : x(x)
  , y(y)
  {};

  IntegerPair(const IntegerPair& other) {
    x = other.x;
    y = other.y;
  }

  inline Bit operator>(const IntegerPair & rhs) const {
    return (x > rhs.x) | ((x == rhs.x) & (y > rhs.y));
  };

  inline IntegerPair operator^=(const IntegerPair& rhs) {
    x ^= rhs.x;
    y ^= rhs.y;
    return (*this);
  };

  inline IntegerPair select(const Bit& sel, const IntegerPair& a) const {
    IntegerPair res(*this);
    res.x = x.select(sel, a.x);
    res.y = y.select(sel, a.y);
    return res;
  }

  inline IntegerPair If(const Bit & sel, const IntegerPair& rhs) const {
    return this->select(sel, rhs);
  }

};
inline IntegerPair If(const Bit & sel, const IntegerPair& a, const IntegerPair& b) {
  IntegerPair res = b;
  return b.If(sel, a);
};

// (a + b) % m, but since a, b already mod, simple check instead
// Somehow sum takes 4 mults, so save when possible
inline Integer AddMod(const Integer& a, const Integer& b, const Integer& m) {
  Integer s = a + b;
  return If(s >= m, s - m, s);
}

struct HeavyEval {
  const int party;

  const CountMin& count_min;
  const HashStorePoly* store;

  const size_t input_range;
  const size_t hash_range;
  const size_t num_hashes;

  // Incremented by 1 since Integers are sometimes signed.
  const size_t input_bits;
  const size_t hash_range_bits;
  // For initial shares, sum modulo big modulus
  const size_t share_bits;
  // (freq) values are max "total count", used just to determine size of values
  const size_t freq_bits;

  Integer* countmin_values = nullptr;

  Integer mod;
  size_t num_values;
  // Sized input_bits
  Integer* values = nullptr;
  // if values is Integer[] from elsewhere, don't double delete
  bool values_is_new = true;
  // Frequencies are Sized freq_bits
  // Tuple <freq, val> for sorting purposes
  IntegerPair* freq_and_vals = nullptr;

  // Total_count only used to determine max size of frequency values
  HeavyEval(const int party, const CountMin& count_min, const size_t total_count)
  : party(party)
  , count_min(count_min)
  , store((HashStorePoly*) count_min.store)
  , input_range(store->get_input_range())
  , hash_range(count_min.cfg.w)
  , num_hashes(count_min.cfg.d)
  , input_bits(store->get_input_bits() + 1)
  , hash_range_bits(LOG2(hash_range) + 1)
  , share_bits(nbits_mod + 1)
  , freq_bits(LOG2(total_count) + 1)
  {
    if (hash_range_bits > input_bits) {
      std::cout << "WARNING: hash range > input range. Count-min not needed. ";
      std::cout << "Behavior will be unreliable" << std::endl;
    }

    mod = Integer(share_bits, fmpz_get_ui(Int_Modulus));
  };

  ~HeavyEval() {
    delete[] countmin_values;
    if (values_is_new) delete[] values;
    delete[] freq_and_vals;
  };

  /*
  Turns shares of count-min to counts as circuit values
  populates countmin_count
  */
  void parse_countmin();

  /*
  Evals i{th} hash on value
  Assume ax + b
  Value: input_bits
  Output: hash_range_bits
  */
  Integer eval_hash(Integer value, const unsigned int i);

  /*
  Gets value from hash_num{th} row of count-min, using index to pick column
  index: hash_range_bits
  output: freq_bits
  */
  Integer get_at_index(const size_t hash_num, Integer index);

  // Array sized num_hashes
  Integer min(Integer* array);

  /* Get frequency estimate of input
  input: input_bits
  output: freq_bits
  */
  Integer get_freq(Integer input);

  // Note: Sorted by freq increasing, so need to take from end of the list.
  void zero_out_dupes();
  void sort_remove_dupes();

  // Populate values with inputs as Integers with input_bits
  // Input is Local secret shares for "party".
  // Num (= num_values) of them.
  void set_values(const fmpz_t* const input_shares, const size_t num);
  // Mostly for testing
  void set_values(const uint64_t* const input_shares, const size_t num);
  // Integer candidate objects from Extract
  void set_values(Integer* inputs, const size_t num);
  // Populate frequencies using values
  void get_frequencies();

  void return_top_K(const size_t K, uint64_t* const topValues,
      uint64_t* const topFreqs);

  void print_params(bool print = true) const {
    if (!print) return;
    std::cout << "HeavyEval params: \n";
    std::cout << "\t Party: " << party << std::endl;
    std::cout << "\t num hashes: " << num_hashes << std::endl;
    std::cout << "\t hash range: " << hash_range << std::endl;
    std::cout << "\t input bits+1: " << input_bits << std::endl;
    std::cout << "\t share bits+1: " << share_bits << std::endl;
    std::cout << "\t freq bits+1: " << freq_bits << std::endl;
    std::cout << "\t hash range bits+1: " << hash_range_bits << std::endl;
  }
  /* DEBUG ONLY */
  void print_countmin(bool print = true) const;
  void print_values(bool print = true) const;
};

/*
Get values (candidates) out of SingleHeavy buckets
*/
struct HeavyExtract {
  const int party;

  const HashStoreBit& store;

  const size_t input_bits;  // for loops, etc

  // for Integer modulus
  // +1 since Integers can be signed
  const size_t value_bits;
  // Initial shares, modulo big modulus
  const size_t share_bits;
  // Base buckets are either x or m-x, so take up full share
  // But signed (and abs) buckets are at most total count
  const size_t freq_bits;

  // Split into B substreams.
  // Each gets D depth SingleHeavy
  // Repeat for Q different splits/hashes
  // Global repeat all R times (same hashes), on different values (halving)
  // R * Q * B * D total buckets
  // Groups of D buckets make value
  // for R * Q * B values
  const size_t R;  // Full copies
  const size_t Q;  // num groups. Num splitting hashes
  const size_t B;  // num substreams
  const size_t D;  // depth

  Integer mod;
  Integer* bucket0;
  Integer* bucket1;
  Bit* cmp;
  Integer* candidates;

  HeavyExtract(const int party, const HashStoreBit& store,
      const size_t R, const size_t Q, const size_t B, const size_t D,
      const size_t total_count)
  : party(party)
  , store(store)
  , input_bits(store.get_input_bits())
  , value_bits(input_bits + 1)
  , share_bits(nbits_mod + 1)
  , freq_bits(LOG2(total_count) + 1)
  , R(R)
  , Q(Q)
  , B(B)
  , D(D)
  {
    mod = Integer(share_bits, fmpz_get_ui(Int_Modulus));
    bucket0 = new Integer[R * Q * B * D];
    bucket1 = new Integer[R * Q * B * D];
    cmp = new Bit[R * Q * B * D];
    candidates = new Integer[R * Q * B];
  };

  ~HeavyExtract() {
    delete[] bucket0;
    delete[] bucket1;
    delete[] cmp;
    delete[] candidates;
  };

  void set_bucket(const int idx, const fmpz_t* const bucket);
  void set_buckets(const fmpz_t* const bucket0, const fmpz_t* const bucket1) {
    set_bucket(0, bucket0);
    set_bucket(1, bucket1);
  }
  // cmp[i] = abs(bucket0[i]) < abs(bucket1[i])
  void bucket_compare();
  // Uses hash inverses on bucket compare values to build candidates
  void extract_candidates();
  Integer* get_candidates() const { return candidates; };

  void print_params(bool print = true) const {
    if (!print) return;
    std::cout << "HeavyExtract params: \n";
    std::cout << "\t input_bits: " << input_bits << "\n";
    std::cout << "\t value_bits (+1): " << value_bits << "\n";
    std::cout << "\t share_bits (+1): " << share_bits << "\n";
    std::cout << "\t freq_bits (+1): " << freq_bits << std::endl;
  }

  /* DEBUG ONLY */
  void print_buckets(bool print = true) const;
  void print_cmp(bool print = true) const;
  void print_candidates(bool print = true) const;
};

/*
NOTE: countmin_shares gets copied into a count-min object
And hence it gets freed at the end of full_heavy_extract.

Standard op:
h = HeavyEval(party, count_min, total_count)
  - Count-min contains secret-shared values
  - total count is # items (sum of frequencies)
Count-min
  parse_countmin()
    - populate countmin_values with de-shared circuit vals
  Can print_countmin past this point
Values
  set_values(int/fmpz* shares, size num_values)
    - Shares of candidate values (potentially with dupes)
    - populates values with actuals with circuit
  can print_values past this point
void get_frequencies()
  - Compute frequency estimates with circuit, populate freq
  - requires both count-min parsed, and values set.
sort_remove_dupes()
  - sort (value/freq) by freq, with dupes zero'd out
return_top_K(size K, uint64_t* topValues, uint64_t* topFreqs)
  - set topValues/Freqs to the top K
  - does circuit reveal at this step
*/

void full_heavy_extract(
    const int server_num, const MultiHeavyConfig cfg,
    const fmpz_t* const bucket0, const fmpz_t* const bucket1,
    flint_rand_t hash_seed_split, flint_rand_t hash_seed_count,
    fmpz_t* const countmin_shares,
    const size_t num_inputs,
    uint64_t* top_values, uint64_t* top_freqs);

#endif
