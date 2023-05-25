#include "eval_heavy.h"

void HeavyEval::parse_countmin() {
  countmin_values = new Integer[num_hashes * hash_range];

  for (unsigned int i = 0; i < num_hashes * hash_range; i++) {
    uint64_t val = fmpz_get_ui(count_min.counts[i]);
    // Can skip this party conditional probably. Here also for clarity
    Integer a(share_bits, party == ALICE ? val : 0, ALICE);
    Integer b(share_bits, party == BOB ? val : 0, BOB);
    countmin_values[i] = (a + b) % mod;
    countmin_values[i].resize(freq_bits, false);
  }
}

Integer HeavyEval::eval_hash(Integer value, const unsigned int i) {
  const size_t mult_bits = input_bits * 2 + 1;

  // Slightly larger to allow for modulo to work right.
  const size_t tmp_size = hash_range_bits + 1;

  Integer a(mult_bits, store->get_coeff(i, 1), PUBLIC);
  Integer b(tmp_size, store->get_coeff(i, 0), PUBLIC);
  Integer hash_range_int(tmp_size, hash_range, PUBLIC);

  // Input range is 2^input_bits, so we do the mod as part of the truncate
  Integer ret = a * value;
  if (input_bits > hash_range_bits) {
    // Should always happen (shrink so count-min gives advantage).
    ret = ret >> (input_bits - hash_range_bits);
  }
  // Quick way to do modulo 2^input_bits
  ret.resize(hash_range_bits-1, false);
  ret.resize(tmp_size, false);
  ret = (ret + b);
  ret = ret % hash_range_int;
  ret.resize(hash_range_bits, false);
  return ret;
}

Integer HeavyEval::get_at_index(const size_t hash_num, Integer index) {
  if (!countmin_values) error_exit("countmin_values not parsed");
  Integer ret(freq_bits, 0, PUBLIC);

  for (unsigned int i = 0; i < hash_range; i++) {
    const size_t offset = hash_num * hash_range + i;

    Integer idx(hash_range_bits, i, PUBLIC);
    ret = If(index == idx, countmin_values[offset], ret);
  }

  // This is secret shared, so need to add/resolve
  // Take in party parameter.
  // Actually, that can be done at the count-min level instead
  // TODO: have two sets of count-min, build a single aggregate one

  return ret;
}

Integer HeavyEval::min(Integer* array) {
  Integer ret = array[0];
  for (unsigned int i = 1; i < num_hashes; i++) {
    ret = If(array[i] < ret, array[i], ret);
  }
  return ret;
}

Integer HeavyEval::get_freq(Integer input) {
  Integer* estimates = new Integer[num_hashes];

  for (unsigned int i = 0; i < num_hashes; i++) {
    Integer hash_val = eval_hash(input, i);
    estimates[i] = get_at_index(i, hash_val);
  }

  Integer ret = min(estimates);
  delete[] estimates;
  return ret;
}

void HeavyEval::zero_out_dupes() {
  if (!frequencies) error_exit("frequencies not set");
  Integer zero_value(input_bits, 0, PUBLIC);
  Integer zero_frequency(freq_bits, 0, PUBLIC);

  /* Can be "more parallel", but not needed.
  but same amount of mult gates, so garbled circuit doesn't care much
  */

  for (unsigned int i = 0; i < num_values - 1; i++) {
    Bit b = (values[i] == values[i+1]);
    values[i] = If(b, zero_value, values[i]);
    frequencies[i] = If(b, zero_frequency, frequencies[i]);
  }
}

void HeavyEval::sort_remove_dupes() {
  if (!frequencies) error_exit("frequencies not set");
  sort(frequencies, num_values, values);
  zero_out_dupes();
  sort(frequencies, num_values, values);
}

void HeavyEval::set_values(Integer* const inputs, const size_t num) {
  num_values = num;
  values = inputs;
  values_is_new = false;
}

void HeavyEval::set_values(const fmpz_t* const input_shares, const size_t num) {
  num_values = num;
  values = new Integer[num_values];
  values_is_new = true;

  for (unsigned int i = 0; i < num_values; i++) {
    uint64_t val = fmpz_get_ui(input_shares[i]);
    // Conditional is overkill, but makes it clearer
    Integer a(share_bits, party == ALICE ? val : 0, ALICE);
    Integer b(share_bits, party == BOB ? val : 0, BOB);
    values[i] = (a + b) % mod;
    values[i].resize(input_bits, false);
  }
}

void HeavyEval::set_values(const uint64_t* const input_shares, const size_t num) {
  num_values = num;
  values = new Integer[num_values];
  values_is_new = true;

  for (unsigned int i = 0; i < num_values; i++) {
    uint64_t val = input_shares[i];
    // Conditional is overkill, but makes it clearer
    Integer a(share_bits, party == ALICE ? val : 0, ALICE);
    Integer b(share_bits, party == BOB ? val : 0, BOB);
    values[i] = (a + b) % mod;
    values[i].resize(input_bits, false);
  }
}

void HeavyEval::get_frequencies() {
  if (!values) error_exit("values not set");
  frequencies = new Integer[num_values];
  for (unsigned int i = 0; i < num_values; i++) {
    frequencies[i] = get_freq(values[i]);
  }
}

void HeavyEval::return_top_K(const size_t K, uint64_t* const topValues, uint64_t* const topFreqs) {
  for (unsigned int i = 0; i < K; i++) {
    topValues[i] = values[num_values - 1 - i].reveal<uint64_t>();
    topFreqs[i] = frequencies[num_values - 1 - i].reveal<uint64_t>();
  }
}

void HeavyEval::print_countmin() const {
  for (unsigned int i = 0; i < num_hashes; i++) {
    for (unsigned int j = 0; j < hash_range; j++) {
      uint64_t x = countmin_values[i * hash_range + j].reveal<uint64_t>();
      if (x)
        std::cout << x << " ";
      else
        std::cout << "_ ";
    }
    std::cout << std::endl;
  }
}

void HeavyEval::print_values() const {
  if (!values) {
    std::cout << "No values" << std::endl;
    return;
  }
  for (unsigned int i = 0; i < num_values; i++) {
    uint64_t v = values[i].reveal<uint64_t>();
    uint64_t f = frequencies ? frequencies[i].reveal<uint64_t>() : 0;
    std::cout << "Value " << i+1 << " / " << num_values << ": " << v;
    if (frequencies) std::cout << ", freq: " << f;
    std::cout << std::endl;
  }
}

void HeavyExtract::set_bucket(const int idx, const fmpz_t* const bucket) {
  Integer* target = idx ? bucket1 : bucket0;

  for (unsigned int i = 0; i < num_buckets; i++) {
    uint64_t val = fmpz_get_ui(bucket[i]);
    Integer a(share_bits, party == ALICE ? val : 0, ALICE);
    Integer b(share_bits, party == BOB ? val : 0, BOB);
    // Since > N/2 is negative, subtract out mod to shift it
    target[i] = (a + b) % mod;
    // No resize, since negatives.
  }
}

void HeavyExtract::bucket_compare() {
  // Absolute value:
  // If x < mod/2, then x (positive)
  // If x > mod/2, then "-x" -> abs(x - mod) = -(x-mod) = mod - x
  Integer mod_half(share_bits, fmpz_get_ui(Int_Modulus)/2, PUBLIC);

  for (unsigned int i = 0; i < num_buckets; i++) {
    Integer abs0 = If(bucket0[i] < mod_half, bucket0[i], mod - bucket0[i]);
    Integer abs1 = If(bucket1[i] < mod_half, bucket1[i], mod - bucket1[i]);
    cmp[i] = abs0 < abs1;
  }
}

void HeavyExtract::extract_candidates() {
  Integer zero(value_bits, 0, PUBLIC);
  Integer two(value_bits, 2, PUBLIC);
  for (unsigned int i = 0; i < num_values; i++) {
    // if (party == ALICE) std::cout << "value idx: " << i << std::endl;
    Integer value(value_bits, 0, PUBLIC);
    for (unsigned int j = 0; j < input_bits; j++) {
      // if (party == ALICE) std::cout << " vec idx: " << j << std::endl;
      Integer vec(value_bits, 0, PUBLIC);
      // const size_t idx = i * input_bits + j;
      for (unsigned int k = 0; k < input_bits; k++) {
        // bool b = cmp[i * input_bits + k].reveal<bool>();
        // if (party == ALICE) std::cout << "   vec[" << k << "] = coeff " << store.get_inv_coeff(i, j, k) << " * cmp " << b << std::endl;
        Integer coeff(value_bits, store.get_inv_coeff(i, j, k), PUBLIC);
        // vec = Multiply Coeff^-1 * bit, which can be just conditional adds
        vec = If(cmp[i * input_bits + k], vec + coeff, vec);
      }
      // Then, use vec mod 2 as ith bit of candidates (times 2^j)
      Integer pow(value_bits, 1ULL<<j, PUBLIC);
      value = If(vec % two == zero, value, value + pow);
    }
    candidates[i] = value;
  }
}

void HeavyExtract::print_buckets() const {
  int64_t m = fmpz_get_ui(Int_Modulus);
  for (unsigned int i = 0; i < num_buckets; i++) {
    int64_t x0 = bucket0[i].reveal<int64_t>();
    int64_t x1 = bucket1[i].reveal<int64_t>();
    x0 -= x0 < m/2 ? 0 : m;
    x1 -= x1 < m/2 ? 0 : m;
    std::cout << "buckets[" << i << "] = " << x0 << " vs " << x1 << std::endl;
  }
}

void HeavyExtract::print_cmp() const {
  int64_t m = fmpz_get_ui(Int_Modulus);
  for (unsigned int i = 0; i < num_buckets; i++) {
    int64_t x0 = bucket0[i].reveal<int64_t>();
    int64_t x1 = bucket1[i].reveal<int64_t>();
    x0 -= x0 < m/2 ? 0 : m;
    x1 -= x1 < m/2 ? 0 : m;
    bool l = cmp[i].reveal<bool>();
    std::cout << "buckets[" << i << "] = |" << x0 << "| ";
    std::cout << (l ? "<":">") << " |" << x1 << "|" << std::endl;
  }
}

void HeavyExtract::print_candidates() const {
  for (unsigned int i = 0; i < num_values; i++) {
    int64_t x = candidates[i].reveal<int64_t>();
    std::cout << "Candidate " << i << " = " << x << std::endl;
  }
}
