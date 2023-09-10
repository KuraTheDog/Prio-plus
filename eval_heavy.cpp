#include "eval_heavy.h"

void HeavyEval::parse_countmin() {
  countmin_values = new Integer[num_hashes * hash_range];

  for (unsigned int i = 0; i < num_hashes * hash_range; i++) {
    uint64_t val = fmpz_get_ui(count_min.counts[i]);
    // Can skip this party conditional probably. Here also for clarity
    Integer a(share_bits, party == ALICE ? val : 0, ALICE);
    Integer b(share_bits, party == BOB ? val : 0, BOB);
    countmin_values[i] = AddMod(a, b, mod);
    countmin_values[i].resize(freq_bits, false);
  }
}

Integer HeavyEval::eval_hash(Integer value, const unsigned int i) {
  const size_t mult_bits = input_bits * 2 + 1;

  // For mod
  // Sizes are +1 to make modulo work better
  // Mult bits has implicit extra +1 that has it work out
  Integer input_range_int(mult_bits, input_range, PUBLIC);
  Integer hash_range_int(hash_range_bits + 1, hash_range, PUBLIC);

  Integer a(mult_bits, store->get_coeff(i, 1), PUBLIC);
  Integer b(hash_range_bits, store->get_coeff(i, 0), PUBLIC);

  // Temp resize value to match for multiply
  value.resize(mult_bits, false);
  Integer ret = a * value;
  value.resize(input_bits, false);

  ret = ret % input_range_int;
  // Grab first bits, from normal hashing
  if (input_bits > hash_range_bits) {
    // Should always happen (shrink so count-min gives advantage).
    ret = ret >> (input_bits - hash_range_bits);
  }
  // Resize to desired, and match b. Above shift stops this from breaking.
  ret.resize(hash_range_bits, false);
  ret = ret + b;
  // Match modulo
  ret.resize(hash_range_bits + 1, false);
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
  Integer* const estimates = new Integer[num_hashes];

  for (unsigned int i = 0; i < num_hashes; i++) {
    Integer hash_val = eval_hash(input, i);
    estimates[i] = get_at_index(i, hash_val);
  }

  Integer ret = min(estimates);
  delete[] estimates;
  return ret;
}

void HeavyEval::zero_out_dupes() {
  if (!freq_and_vals) error_exit("frequencies not set");
  IntegerPair zero_pair(Integer(freq_bits, 0, PUBLIC), Integer(input_bits, 0, PUBLIC));
  Integer zero_value(input_bits, 0, PUBLIC);
  Integer zero_frequency(freq_bits, 0, PUBLIC);

  /* Can be "more parallel", but not needed.
  but same amount of mult gates, so garbled circuit doesn't care much
  */

  for (unsigned int i = 0; i < num_values - 1; i++) {
    // Since values have unique freq, just have to check if values are the same
    Bit b = (freq_and_vals[i].y == freq_and_vals[i+1].y);
    freq_and_vals[i] = If(b, zero_pair, freq_and_vals[i]);
  }
}

void HeavyEval::sort_remove_dupes() {
  if (!freq_and_vals) error_exit("frequencies not set");
  sort(freq_and_vals, num_values);
  // std::cout << "sorted values: " << std::endl; print_values();
  zero_out_dupes();
  // std::cout << "sorted values, dupes zeroed out: " << std::endl; print_values();
  sort(freq_and_vals, num_values);
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
    values[i] = AddMod(a, b, mod);
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
    values[i] = AddMod(a, b, mod);
    values[i].resize(input_bits, false);
  }
}

void HeavyEval::get_frequencies() {
  if (!values) error_exit("values not set");
  freq_and_vals = new IntegerPair[num_values];
  for (unsigned int i = 0; i < num_values; i++) {
    freq_and_vals[i].y = values[i];
    freq_and_vals[i].x = get_freq(values[i]);
  }
}

void HeavyEval::return_top_K(const size_t K,
    uint64_t* const topValues, uint64_t* const topFreqs) {
  for (unsigned int i = 0; i < K; i++) {
    IntegerPair pair = freq_and_vals[num_values - 1 - i];
    topValues[i] = pair.y.reveal<uint64_t>();
    topFreqs[i] = pair.x.reveal<uint64_t>();
  }
}

void HeavyEval::print_countmin(bool print) const {
  for (unsigned int i = 0; i < num_hashes; i++) {
    for (unsigned int j = 0; j < hash_range; j++) {
      uint64_t x = countmin_values[i * hash_range + j].reveal<uint64_t>();
      if (print) {
        if (x)
          std::cout << x << " ";
        else
          std::cout << "_ ";
      }
    }
    if (print) std::cout << std::endl;
  }
}

void HeavyEval::print_values(bool print) const {
  if (!values) {
    std::cout << "No values" << std::endl;
    return;
  }
  for (unsigned int i = 0; i < num_values; i++) {
    uint64_t v, f;
    if (freq_and_vals) {
      v = freq_and_vals[i].y.reveal<uint64_t>();
      f = freq_and_vals[i].x.reveal<uint64_t>();
    } else {
      v = values[i].reveal<uint64_t>();
    }
    if (print) {
      std::cout << "Value " << i+1 << " / " << num_values << ": " << v;
      if (freq_and_vals) std::cout << ", freq: " << f;
      std::cout << std::endl;
    }
  }
}

// TODO: Somehow this is the slowest part? Just a lot of declares
void HeavyExtract::set_bucket(const int idx, const fmpz_t* const bucket) {
  Integer* const target = idx ? bucket1 : bucket0;

  for (unsigned int i = 0; i < R * Q * B * D; i++) {
    uint64_t val = fmpz_get_ui(bucket[i]);
    Integer a(share_bits, party == ALICE ? val : 0, ALICE);
    Integer b(share_bits, party == BOB ? val : 0, BOB);
    target[i] = AddMod(a, b, mod);
  }
}

void HeavyExtract::bucket_compare() {
  // Absolute value:
  // If x < mod/2, then x (positive)
  // If x > mod/2, then "-x" -> abs(x - mod) = -(x-mod) = mod - x
  Integer mod_half(share_bits, fmpz_get_ui(Int_Modulus)/2, PUBLIC);

  for (unsigned int i = 0; i < R * Q * B * D; i++) {
    Integer abs0 = If(bucket0[i] < mod_half, bucket0[i], mod - bucket0[i]);
    Integer abs1 = If(bucket1[i] < mod_half, bucket1[i], mod - bucket1[i]);
    abs0.resize(freq_bits, false);
    abs1.resize(freq_bits, false);
    cmp[i] = abs0 < abs1;
  }
}

void HeavyExtract::extract_candidates() {
  Integer zero(value_bits, 0, PUBLIC);
  Integer one(value_bits, 1, PUBLIC);
  Integer two(value_bits, 2, PUBLIC);

  unsigned int val_idx = 0;
  // Global repeats R
  for (unsigned int r = 0; r < R; r++) {
    // Row repeat Q. Determines hash
    for (unsigned int q = 0; q < Q; q++) {
      // Which substream. Builds value for R*Q*B
      for (unsigned int b = 0; b < B; b++) {
        Integer value(value_bits, 0, PUBLIC);
        // Solving over bits of value
        for (unsigned int bit = 0; bit < input_bits; bit++) {
          // Since do value % 2 at the end, just short circuit
          Bit b = Bit(0, PUBLIC);
          // Over depth, inverting and adding (mod 2)
          for (unsigned int d = 0; d < D; d++) {
            const unsigned int cmp_idx = val_idx * D + d;
            // Coeffs are 0 or 1, and fixed
            // So we can reduce circuit size by conditional on public coeff
            if (store.get_inv_coeff(q, bit, d) == 1) {
              b = If(cmp[cmp_idx], !b, b);
            }
          }
          // Can store powers for reuse, but probably fine since public const
          Integer pow(value_bits, 1ULL << bit, PUBLIC);
          value = If(b, value + pow, value);
        }
        candidates[val_idx] = value;
        val_idx++;
      }
    }
  }
}

void HeavyExtract::print_buckets(bool print) const {
  int64_t m = fmpz_get_ui(Int_Modulus);
  for (unsigned int i = 0; i < R * Q * B * D; i++) {
    int64_t x0 = bucket0[i].reveal<int64_t>();
    int64_t x1 = bucket1[i].reveal<int64_t>();
    x0 -= x0 < m/2 ? 0 : m;
    x1 -= x1 < m/2 ? 0 : m;
    if (print)
      std::cout << "buckets[" << i << "] = " << x0 << " vs " << x1 << std::endl;
  }
}

void HeavyExtract::print_cmp(bool print) const {
  int64_t m = fmpz_get_ui(Int_Modulus);
  for (unsigned int i = 0; i < R * Q * B * D; i++) {
    int64_t x0 = bucket0[i].reveal<int64_t>();
    int64_t x1 = bucket1[i].reveal<int64_t>();
    x0 -= x0 < m/2 ? 0 : m;
    x1 -= x1 < m/2 ? 0 : m;
    bool l = cmp[i].reveal<bool>();
    if (print) {
      std::cout << "buckets[" << i << "] = |" << x0 << "| ";
      std::cout << (l ? "<":">") << " |" << x1 << "|" << std::endl;
    }
  }
}

void HeavyExtract::print_candidates(bool print) const {
  for (unsigned int i = 0; i < R * Q * B; i++) {
    int64_t x = candidates[i].reveal<int64_t>();
    if (print)
      std::cout << "Candidate " << i << " / " << R * Q * B << " = " << x << std::endl;
  }
}

void full_heavy_extract(
    const int server_num, const MultiHeavyConfig cfg,
    const fmpz_t* const bucket0, const fmpz_t* const bucket1,
    flint_rand_t hash_seed_split, flint_rand_t hash_seed_count,
    fmpz_t* const countmin_shares,
    const size_t num_inputs,
    uint64_t* top_values, uint64_t* top_freqs) {
  const int party = server_num + 1;

  HashStoreBit hash_split(cfg.Q, cfg.D, cfg.num_bits, 2, hash_seed_split);

  HeavyExtract ex(party, hash_split, cfg.R, cfg.Q, cfg.B, cfg.D, num_inputs);
  // ex.print_params();

  ex.set_buckets(bucket0, bucket1);
  // std::cout << "buckets: " << std::endl; ex.print_buckets();

  ex.bucket_compare();
  // std::cout << "compare: " << std::endl; ex.print_cmp();

  ex.extract_candidates();
  // std::cout << "candidates: " << std::endl; ex.print_candidates();

  CountMin count_min(cfg.countmin_cfg);
  count_min.setStore(cfg.num_bits, hash_seed_count);
  count_min.counts = countmin_shares;
  // count_min.print();

  HeavyEval ev(party, count_min, num_inputs);
  // ev.print_params();

  ev.parse_countmin();
  // std::cout << "countmin: " << std::endl; ev.print_countmin();

  ev.set_values(ex.get_candidates(), cfg.R * cfg.Q * cfg.B);
  // std::cout << "set values: " << std::endl; ev.print_values();

  ev.get_frequencies();
  // std::cout << "values and freqs: " << std::endl; ev.print_values();
  ev.sort_remove_dupes();
  // std::cout << "values and freqs, sort deduped: " << std::endl; ev.print_values();
  ev.return_top_K(cfg.K, top_values, top_freqs);
}
