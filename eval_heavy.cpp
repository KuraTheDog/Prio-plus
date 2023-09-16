#include <set>

#include "eval_heavy.h"
#include "net_share.h"

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

void HeavyEval::get_hashes() {
  if (!values) error_exit("values not set");
  hashes = new Integer[num_values * num_hashes];

  for (unsigned int i = 0; i < num_values; i++)
    for (unsigned int j = 0; j < num_hashes; j++)
      hashes[i * num_hashes + j] = eval_hash(values[i], j);
}

uint64_t HeavyEval::get_freq(Integer input) {
  Integer* const estimates = new Integer[num_hashes];

  for (unsigned int i = 0; i < num_hashes; i++) {
    Integer hash_val = eval_hash(input, i);
    estimates[i] = get_at_index(i, hash_val);
  }

  Integer ret = min(estimates);
  delete[] estimates;
  return ret.reveal<uint64_t>();
}

void HeavyEval::get_frequencies() {
  if (!values) error_exit("values not set");
  if (!hashes) error_exit("hashes not set");
  freq_and_vals = new FreqVal[num_values];

  Integer* const estimates = new Integer[num_hashes];

  for (unsigned int i = 0; i < num_values; i++) {
    for (unsigned int j = 0; j < num_hashes; j++) {
      estimates[j] = get_at_index(j, hashes[i * num_hashes + j]);
    }
    Integer freq = min(estimates);
    freq_and_vals[i] = FreqVal(freq.reveal<uint64_t>(), values[i]);
  }

  delete[] estimates;

  delete[] hashes; hashes = nullptr;
  delete[] countmin_values; countmin_values = nullptr;
}

void HeavyEval::sort_remove_dupes(const size_t K) {
  topK = new pair[K];

  // Sort via freq, ignore values
  std::sort(freq_and_vals, freq_and_vals + num_values,
      [](FreqVal a, FreqVal b) { return a.first > b.first; });

  size_t offset = 0;
  size_t unique = 0;
  while (unique < K) {
    // Get curr
    FreqVal curr = freq_and_vals[offset];
    uint64_t curr_freq = curr.first;
    /* Edge case: Not enough candidates have non-zero freq. Just stop early
    Doesn't leak to reveal extra values that have 0
    - i.e. anything with 0 freq is random values, so no leak
    - and it would be revealed that some < K are all anyways
    - so don't need to "stop top K early" in final return
    */
    if (curr_freq == 0) break;
    // Get all with same freq as curr
    size_t count = 0;
    while (curr.first == curr_freq) {
      count += 1;
      curr = freq_and_vals[offset + count];
    }

    // Reveal values with the freq. No leakage
    std::set<uint64_t> vals;
    for (unsigned int i = offset; i < offset + count; i++) {
      vals.insert(freq_and_vals[i].second.reveal<uint64_t>());
    }

    // Add unique values to top K
    size_t idx = unique;
    for (auto v : vals) {
      topK[idx] = {v, curr_freq};
      idx++;
      if (idx == K) break;  // Found enough, stop
    }

    unique += vals.size();
    offset += count;
  }
}

void HeavyEval::set_values(Integer* const inputs, const size_t num) {
  num_values = num;
  values = inputs;
  values_is_new = false;
}

void HeavyEval::set_values(const fmpz_t* const shares, const size_t num) {
  num_values = num;
  values = new Integer[num_values];
  values_is_new = true;

  for (unsigned int i = 0; i < num_values; i++) {
    uint64_t val = fmpz_get_ui(shares[i]);
    // Conditional is overkill, but makes it clearer
    Integer a(share_bits, party == ALICE ? val : 0, ALICE);
    Integer b(share_bits, party == BOB ? val : 0, BOB);
    values[i] = AddMod(a, b, mod);
    values[i].resize(input_bits, false);
  }
}

void HeavyEval::set_values(const uint64_t* const shares, const size_t num) {
  num_values = num;
  values = new Integer[num_values];
  values_is_new = true;

  for (unsigned int i = 0; i < num_values; i++) {
    uint64_t val = shares[i];
    // Conditional is overkill, but makes it clearer
    Integer a(share_bits, party == ALICE ? val : 0, ALICE);
    Integer b(share_bits, party == BOB ? val : 0, BOB);
    values[i] = AddMod(a, b, mod);
    values[i].resize(input_bits, false);
  }
}

void HeavyEval::set_values_clear(const fmpz_t* const shares, const size_t num,
    const int who) {
  num_values = num;
  values = new Integer[num_values];
  values_is_new = true;

  for (unsigned int i = 0; i < num_values; i++)
    values[i] = Integer(input_bits, fmpz_get_ui(shares[i]), who);
}

void HeavyEval::set_values_clear(const uint64_t* const shares, const size_t num,
    const int who) {
  num_values = num;
  values = new Integer[num_values];
  values_is_new = true;

  for (unsigned int i = 0; i < num_values; i++)
    values[i] = Integer(input_bits, shares[i], who);
}

void HeavyEval::set_hashes(const fmpz_t* const shares) {
  hashes = new Integer[num_values * num_hashes];

  for (unsigned int i = 0; i < num_values * num_hashes; i++) {
    uint64_t val = fmpz_get_ui(shares[i]);
    // Conditional is overkill, but makes it clearer
    Integer a(share_bits, party == ALICE ? val : 0, ALICE);
    Integer b(share_bits, party == BOB ? val : 0, BOB);
    hashes[i] = AddMod(a, b, mod);
    hashes[i].resize(hash_range_bits, false);
  }
}

void HeavyEval::set_hashes(const uint64_t* const shares) {
  hashes = new Integer[num_values * num_hashes];

  for (unsigned int i = 0; i < num_values * num_hashes; i++) {
    uint64_t val = shares[i];
    // Conditional is overkill, but makes it clearer
    Integer a(share_bits, party == ALICE ? val : 0, ALICE);
    Integer b(share_bits, party == BOB ? val : 0, BOB);
    hashes[i] = AddMod(a, b, mod);
    hashes[i].resize(hash_range_bits, false);
  }
}

void HeavyEval::set_hashes_clear(const fmpz_t* const shares, const int who) {
  hashes = new Integer[num_values * num_hashes];
  for (unsigned int i = 0; i < num_values * num_hashes; i++)
    hashes[i] = Integer(hash_range_bits, fmpz_get_ui(shares[i]), who);
}

void HeavyEval::set_hashes_clear(const uint64_t* const shares, const int who) {
  hashes = new Integer[num_values * num_hashes];
  for (unsigned int i = 0; i < num_values * num_hashes; i++)
    hashes[i] = Integer(hash_range_bits, shares[i], who);
}

void HeavyEval::return_top_K(const size_t K,
    uint64_t* const topValues, uint64_t* const topFreqs) {
  if (!topK) error_exit("top K not found yet");
  for (unsigned int i = 0; i < K; i++) {
    std::tie(topValues[i], topFreqs[i]) = topK[i];
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
      v = freq_and_vals[i].second.reveal<uint64_t>();
      f = freq_and_vals[i].first;
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

  delete[] bucket0; bucket0 = nullptr;
  delete[] bucket1; bucket1 = nullptr;
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
          Bit parity = Bit(0, PUBLIC);
          // Over depth, inverting and adding (mod 2)
          for (unsigned int d = 0; d < D; d++) {
            const unsigned int cmp_idx = val_idx * D + d;
            // Coeffs are 0 or 1, and fixed
            // So we can reduce circuit size by conditional on public coeff
            if (store.get_inv_coeff(q, bit, d) == 1) {
              parity = If(cmp[cmp_idx], !parity, parity);
            }
          }
          // Can store powers for reuse, but probably fine since public const
          Integer pow(value_bits, 1ULL << bit, PUBLIC);
          value = If(parity, value + pow, value);
        }
        candidates[val_idx] = value;
        val_idx++;
      }
    }
  }

  delete[] cmp; cmp = nullptr;
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

bool* bucket_compare_clear(
    const int serverfd, const size_t N,
    const fmpz_t* const bucket0, const fmpz_t* const bucket1) {
  // const size_t N = R * Q * B * D;
  bool* const cmp = new bool[N];
  // Reveal buckets. Double send instead of extra copy into buff
  fmpz_t* buff0; new_fmpz_array(&buff0, N);
  fmpz_t* buff1; new_fmpz_array(&buff1, N);
  std::thread t_send([serverfd, bucket0, bucket1, N]() {
    send_fmpz_batch(serverfd, bucket0, N);
    send_fmpz_batch(serverfd, bucket1, N);
  });
  std::thread t_recv([serverfd, &buff0, &buff1, N]() {
    recv_fmpz_batch(serverfd, buff0, N);
    recv_fmpz_batch(serverfd, buff1, N);
  });
  t_send.join();
  t_recv.join();

  // Not using util, to reuse half declare
  fmpz_t half; fmpz_init(half); fmpz_cdiv_q_ui(half, Int_Modulus, 2);
  for (unsigned int i = 0; i < N; i++) {
    fmpz_mod_add(buff0[i], bucket0[i], buff0[i], mod_ctx);
    fmpz_mod_add(buff1[i], bucket1[i], buff1[i], mod_ctx);
    // Abs: if > M/2, then x = x-M < 0, so |x| = M - x
    if (fmpz_cmp(buff0[i], half) > 0) fmpz_sub(buff0[i], Int_Modulus, buff0[i]);
    if (fmpz_cmp(buff1[i], half) > 0) fmpz_sub(buff1[i], Int_Modulus, buff1[i]);
    cmp[i] = (fmpz_cmp(buff0[i], buff1[i]) < 0);
  }
  clear_fmpz_array(buff0, N);
  clear_fmpz_array(buff1, N);
  fmpz_clear(half);

  return cmp;
}

void extract_candidates_clear(
    const size_t R, const size_t Q, const size_t B, const size_t D,
    const size_t input_bits, const HashStoreBit& store,
    const bool* const cmp, uint64_t* const candidates) {

  // Compute candidates, size R * Q * B
  unsigned int val_idx = 0;
  // Global repeats R
  for (unsigned int r = 0; r < R; r++) {
    // Row repeat Q. Determines hash
    for (unsigned int q = 0; q < Q; q++) {
      // Which substream. Builds value for R*Q*B
      for (unsigned int b = 0; b < B; b++) {
        uint64_t value = 0;
        // Solving over bits of value
        for (unsigned int bit = 0; bit < input_bits; bit++) {
          // Since do value % 2 at the end, just short circuit
          bool parity = false;
          // Over depth, inverting and adding (mod 2)
          for (unsigned int d = 0; d < D; d++) {
            const unsigned int cmp_idx = val_idx * D + d;
            // Coeffs are 0 or 1, and fixed
            if (cmp[cmp_idx] and store.get_inv_coeff(q, bit, d) == 1) {
              parity ^= 1;
            }
          }
          // += parity << bit
          if (parity)
            value += (1ULL << bit);
        }
        candidates[val_idx] = value;
        val_idx++;
      }
    }
  }
}

void top_k_extract_garbled(
    const int server_num, const MultiHeavyConfig cfg,
    const fmpz_t* const bucket0, const fmpz_t* const bucket1,
    flint_rand_t hash_seed_split, flint_rand_t hash_seed_count,
    fmpz_t* const countmin_shares, const size_t num_inputs,
    uint64_t* const top_values, uint64_t* const top_freqs) {
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

  ev.get_hashes();

  ev.get_frequencies();
  // std::cout << "values and freqs: " << std::endl; ev.print_values();

  ev.sort_remove_dupes(cfg.K);

  ev.return_top_K(cfg.K, top_values, top_freqs);
}

void top_k_extract_mixed(
    const int server_num, const int serverfd, const MultiHeavyConfig cfg,
    const fmpz_t* const bucket0, const fmpz_t* const bucket1,
    flint_rand_t hash_seed_split, flint_rand_t hash_seed_count,
    fmpz_t* const countmin_shares, const size_t num_inputs,
    uint64_t* const top_values, uint64_t* const top_freqs) {
  const int party = server_num + 1;
  /* Currently reveals all candidates in the clear
      - Leaks possibly popular inputs
      - could also include random, or not popular. Freq still hidden 
  TODO: abs_cmp get bit, eval hash on it.
  */
  const size_t N = cfg.R * cfg.Q * cfg.B * cfg.D;
  const size_t total_candidates = cfg.R * cfg.Q * cfg.B;
  const size_t input_bits = cfg.num_bits;

  CountMin count_min(cfg.countmin_cfg);
  count_min.setStore(cfg.num_bits, hash_seed_count);
  count_min.counts = countmin_shares;
  const size_t num_hashes = count_min.cfg.d;

  HashStoreBit hash_split(cfg.Q, cfg.D, cfg.num_bits, 2, hash_seed_split);

  // Candidates and bucket compares in clear. 
  bool* const cmp = bucket_compare_clear(serverfd, N, bucket0, bucket1);

  uint64_t* const candidates = new uint64_t[total_candidates];
  extract_candidates_clear(cfg.R, cfg.Q, cfg.B, cfg.D, input_bits,
      hash_split, cmp, candidates);
  delete[] cmp;

  // Since candidates in the clear, can de-dupe them now.
  // Means less garbled querying
  std::set<uint64_t> unique_candidates;
  for (unsigned int i = 0; i < total_candidates; i++) {
    unique_candidates.insert(candidates[i]);
  }
  const size_t num_unique_candidates = unique_candidates.size();
  std::cout << "unique candidates: " << unique_candidates.size();
  std::cout << " / " << total_candidates << std::endl;

  // Since in clear, only set for server 0 to get shares
  if (server_num == 0) {
    size_t i = 0;
    for (auto v : unique_candidates) {
      candidates[i] = v;
      i++;
    }
  } else {
    memset(candidates, 0, total_candidates * sizeof(uint64_t));
  }

  // Since in clear, only do for server 0 to get shares
  fmpz_t* hashes; new_fmpz_array(&hashes, num_unique_candidates * num_hashes);
  if (server_num == 0) {
    for (unsigned int i = 0; i < num_unique_candidates; i++) {
      for (unsigned int j = 0; j < num_hashes; j++) {
        count_min.store->eval(j, candidates[i], hashes[i * num_hashes + j]);
      }
    }
  }


  HeavyEval ev(party, count_min, num_inputs);
  // ev.print_params();

  ev.parse_countmin();
  // std::cout << "countmin: " << std::endl; ev.print_countmin();

  ev.set_values_clear(candidates, num_unique_candidates);
  delete[] candidates;
  // std::cout << "set values: " << std::endl; ev.print_values();

  ev.set_hashes_clear(hashes);
  clear_fmpz_array(hashes, num_unique_candidates * num_hashes);

  ev.get_frequencies();
  // std::cout << "values and freqs: " << std::endl; ev.print_values();

  ev.sort_remove_dupes(cfg.K);

  ev.return_top_K(cfg.K, top_values, top_freqs);
}
