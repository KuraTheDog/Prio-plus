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

Note: Seem to need to set up sockets beforehand.
*/

/*
Pieces can be abstracted out
But for now, in one struct so easier to access various attributes.
Mostly one-time, but can call with multiple K I guess.
Can theoretically be one big function with sub-functions.
Currently split for clarity of use / debugging, and for replacements

Note: emp::ALICE is 1, emp::BOB is 2. (Using server_num + 1)
*/

/* Some mult gate costs, for n bit numbers:

Delares: free
-a: n - 1
a +- b: n - 1 (since sign bit)
a * b: n^2
a / b: n^2 + 6n - 4
a % b: n^2 + 6n - 4
a < b: n
a == b : n
If(bit, a, b): n
>>, << const: free
a <<, >> b: 4 n

Bits:
!, ^: free
&, |: 1
If(, bit, bit): 1

Notes:
- Always rounds up to next mult of 4
- No optimization, re=computes. So storing intermediates better if reused
- Doesn't matter if values are public or individual
*/

typedef std::pair<uint64_t, Integer> FreqVal;

// (a + b) % m, but since a, b already mod, simple check instead
// Gate cost: 2 * (nbits-1) + 2 * nbits = 4 * nbits - 2
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
  Integer* hashes = nullptr;
  FreqVal* freq_and_vals = nullptr;
  typedef std::pair<uint64_t, uint64_t> pair;
  pair* topK = nullptr;

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
    delete[] hashes;
    delete[] freq_and_vals;
    delete[] topK;
  };

  /*
  Turns shares of count-min to counts as circuit values
  populates countmin_count
  Gates: 4 * (share_bits - 1)
  */
  void parse_countmin();

  /*
  Evals i{th} hash on value
  Assume ax + b
  Value: input_bits
  Output: hash_range_bits
  Gates: (mul, mod, add, mod)
    - m = mult bits = 2 * input_bits + 1, h = hash_range
    - m^2 + (m^2 + 6m - 4) + h-1 + h^2 + 8h + 3
    = 2m^2 + 6m + h^2 + 9h - 2
    = 8 inp^2 + 20 inp + hash^2 + 9 hash + 6
  */
  Integer eval_hash(Integer value, const unsigned int i);

  /*
  Gets value from hash_num{th} row of count-min, using index to pick column
  index: hash_range_bits
  output: freq_bits
  Gates: 2 * freq_bits * hash_range
  */
  Integer get_at_index(const size_t hash_num, Integer index);

  // Array sized num_hashes
  // Gates: 2 * nbits * num_hashes
  Integer min(Integer* array);

  // Populate hashes
  // Gates: num_values * num_hashes * eval_hash
  void get_hashes();

  /* Get frequency estimate of input from input.
  Unused, replaced by populating hashes
  input: input_bits
  output: in the clear, since value still hidden
  Gates: num_hashes * (eval + get_at) + min
  */
  uint64_t get_freq(Integer input);

  /* Populate frequencies using hashes
  Also clears out hashes and count-min
  Gates: num_values * (num_hashes * get_at + min)
  */
  void get_frequencies();

  /* Reveals just frequencies, but not values, which is minimal leakage
     No gates used, since ops on revealed only
     De-dupes via revealing valuse at freq.
     Stops when >= K unique values found.
     Also "leaks" values that tie Kth largest, which is fine.
  */
  void sort_remove_dupes(const size_t K);

  // Populate values with inputs as Integers with input_bits
  // Input is Local secret shares for "party".
  // Gates (non-int): (4 * share_bits - 2) * num
  // Sets num_values = Num
  void set_values(const fmpz_t* const shares, const size_t num);
  void set_values(const uint64_t* const shares, const size_t num);
  // Integer candidate objects from Extract
  void set_values(Integer* inputs, const size_t num);
  // If not shares, don't need mod. Assumes "who" has the values, other has 0
  void set_values_clear(const fmpz_t* const shares, const size_t num,
      const int who = ALICE);
  void set_values_clear(const uint64_t* const shares, const size_t num,
      const int who = ALICE);

  // Populate hashes
  // Gates: num_values * num_hashes * AddMod
  void set_hashes(const fmpz_t* const shares);
  void set_hashes(const uint64_t* const shares);
  // If not shares, don't need mod. Assumes "who" has the values, other has 0
  void set_hashes_clear(const fmpz_t* const shares, const int who = ALICE);
  void set_hashes_clear(const uint64_t* const shares, const int who = ALICE);

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

  // Gates: Shares, so (4 * share_bits - 2) * RQBD
  void set_bucket(const int idx, const fmpz_t* const bucket);
  void set_buckets(const fmpz_t* const bucket0, const fmpz_t* const bucket1) {
    set_bucket(0, bucket0);
    set_bucket(1, bucket1);
  }
  // cmp[i] = abs(bucket0[i]) < abs(bucket1[i])
  // Gates: (6 * share_bits + freq_bits - 2) * RQBD
  void bucket_compare();
  /* Uses hash inverses on bucket compare values to build candidates
    Gates:
      - Scales with # 1s in inverse matrix N
      - RQB(input_bits) (D * N + value_bits)
  */
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
Possible things to toggle.
Note that the slowest thing by far is just setting up garble for bucket

Buckets:
- clear: reveal bucket values
- precompute: using abs_cmp. Gives arithmetic shares of cmp
  - TODO: figure out if can do boolean shares
- Garble

Extract:
- clear: clear cmp in, clear value out
- secure: [cmp] in, [value] out
  - TODO: Requires cmp ^ parity
  - Seems straightforward if [cmp] is bool. TODO: see if LSB can return bool
  - if [cmp]^A, then need mult round for xor
    - xor(list), so multiple rounds, oof
  - Final is sum_{parity_i * 1<<i} essentially
  - so b2A round to get final values
- garble

Count-min
- full clear: reveal struct
- Garble: keep struct hidden. Useful

Final sort
- clear
- garble
- freq in clear
- no real other way to do efficiently
*/
// Full clear cmp, also reveals bucket value
bool* bucket_compare_clear(const int serverfd, const size_t N,
    const fmpz_t* const bucket0, const fmpz_t* const bucket1);

// Abs_cmp, in correlated
// TODO: Investigate if can have it return bool efficiently
//   Doable, but cmp_bit is a lot of rounds, and can't be compacted for bool

// Given clear cmp, eval candidates
void extract_candidates_clear(
    const size_t R, const size_t Q, const size_t B, const size_t D,
    const size_t input_bits, const HashStoreBit& store,
    const bool* const cmp, uint64_t* const candidates);

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

void top_k_extract_garbled(
    const int server_num, const MultiHeavyConfig cfg,
    const fmpz_t* const bucket0, const fmpz_t* const bucket1,
    flint_rand_t hash_seed_split, flint_rand_t hash_seed_count,
    fmpz_t* const countmin_shares, const size_t num_inputs,
    uint64_t* const top_values, uint64_t* const top_freqs);

#endif
