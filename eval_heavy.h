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

TODO: bucket evaluation with cmp too
*/

/*
Pieces can be abstracted out
But for now, in one struct so easier to access various attributes.
Mostly one-time, but can call with multiple K I guess.
Can theoretically be one big function with sub-functions.
Currently split for clarity of use / debugging.

There may be better ways of de-duping, but this should be fine.
 - challenge is to keep it oblivious. same work if all dupes and no dupes.
 - Can't track how many distinct values seen
 - can't pre-prune list of candidates.

Note: emp::ALICE is 1, emp::BOB is 2. Not quite server_num.
*/

/*
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

struct HeavyEval {
  const int party;

  const CountMin& count_min;
  const HashStorePoly* store;

  const size_t hash_range;
  const size_t num_hashes;
  const size_t total_count;

  // Incremented by 1 since Integers are sometimes signed.
  const size_t input_bits;
  const size_t hash_range_bits;
  // For initial shares, sum modulo big modulus
  const size_t share_bits;
  // (freq) values are max "total count"
  const size_t freq_bits;

  Integer* countmin_values = nullptr;

  size_t num_values;
  // Sized input_bits
  Integer* values = nullptr;
  // Sized freq_bits
  Integer* frequencies = nullptr;

  HeavyEval(const int party, CountMin& count_min, const size_t total_count)
  : party(party)
  , count_min(count_min)
  , store((HashStorePoly*) count_min.store)
  , hash_range(count_min.cfg.w)
  , num_hashes(count_min.cfg.d)
  , total_count(total_count)
  , input_bits(store->get_input_bits() + 1)
  , hash_range_bits(LOG2(hash_range) + 1)
  , share_bits(nbits_mod + 1)
  , freq_bits(LOG2(total_count) + 1)
  {
    std::cout << "HeavyEval params: \n";
    std::cout << "\t Party: " << party << std::endl;
    std::cout << "\t num hashes: " << num_hashes << std::endl;
    std::cout << "\t hash range: " << hash_range << std::endl;
    std::cout << "\t total_count: " << total_count << std::endl;
    std::cout << "\t input bits+1: " << input_bits << std::endl;
    std::cout << "\t share bits+1: " << share_bits << std::endl;
    std::cout << "\t freq bits+1: " << freq_bits << std::endl;
    std::cout << "\t hash range bits+1: " << hash_range_bits << std::endl;

    if (hash_range_bits > input_bits) {
      std::cout << "WARNING: hash range > input range. Count-min not needed. ";
      std::cout << "Behavior will be unreliable" << std::endl;
    }
  };

  ~HeavyEval() {
    delete[] countmin_values;
    delete[] values;
    delete[] frequencies;
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

  void return_top_K(const size_t K, uint64_t* const topValues, uint64_t* const topFreqs);

  /* DEBUG ONLY */
  void print_countmin();
  void print_values();
};

#endif
