#undef NDEBUG
#include <iostream>

#include "utils_test_connect.h"
#include "../constants.h"
#include "../eval_heavy.h"
#include "../hash.h"
#include "../heavy.h"
#include "../ot.h"

#define SERVER0_IP "127.0.0.1"
#define SERVER1_IP "127.0.0.1"

/*
Integer(nbits, val, party):

reveal<type>(who gets it)
- Pub reveal: both have
- Alice reveal: Alice has
- Bob reveal: bob has
- Independent of original value owner

All math is mod nbits
Want nbits >> modulus.

Need to kinda just do secret sharing on it's own.
[v] split into Va party Alice, Vb party Bob
and v = (Va + Vb) % mod

Build "true" values as [(Va + Vb) % mod], and do math on those.

There's no easy "share mult", just write it out.

Constants: make public value.
*/

void test_compare(const int party) {

  Integer a(4, 3, ALICE);
  Integer b(4, 2, BOB);
  Bit res = a > b;

  std::cout << "Party " << party << ", ALICE larger?\t";
  std::cout << res.reveal<bool>() << std::endl;
}

void test_mod(const int party) {
  Integer a(4, 3, ALICE);
  Integer b(4, 2, BOB);
  Integer m(4, 4, PUBLIC);

  // V1: 40 mult gates
  // Old full: 449740
  // Integer res = (a + b) % m;

  /*V2:
  16: s, if(s>m, s-m,s)
  16: s, 0, s - if(s>m, m, 0)
  20: s, if(s>m, s-m, a+b)
  24: if(a+b>m, a+b, 0)
  Somehow add takes 4 (per input). So save when possible
    Re-use of values helps
    since bitwise representation, gets costly

  a + b - m itself is 8 already, half as much
  New full: 65780, nearly 7x smaller
  */
  // Integer s = a + b;
  // Integer res = If(s > m, s-m, s);
  // Integer res = a + b - m;
  Integer res = AddMod(a, b, m);

  std::cout << "3 + 4 = " << res.reveal<uint64_t>() << " mod 4\n";

}

void test_hash(const int party, flint_rand_t hash_seed) {

  CountMinConfig cfg(.1, .4);  // 3 hashes range 7
  const size_t input_bits = 4;  // At least >= range bits
  CountMin count_min(cfg);
  count_min.init();
  count_min.setStore(input_bits, hash_seed);
  HeavyEval heavy_eval(party, count_min, 10);

  size_t input_size = 1ULL << input_bits;

  fmpz_t out; fmpz_init(out);

  for (unsigned int i = 0; i < input_size; i++) {
    Integer x(input_bits + 1, i, PUBLIC);
    for (unsigned int j = 0; j < heavy_eval.num_hashes; j++) {
      Integer ret = heavy_eval.eval_hash(x, j);
      int64_t garble_out = ret.reveal<int64_t>();
      heavy_eval.count_min.store->eval(j, i, out);
      int64_t hash_out = fmpz_get_ui(out);
      // if (garble_out != hash_out) {
      //   std::cout << "failed on hash_" << j << "(" << i << ")" << std::endl;
      //   count_min.store->print_hash(j);
      //   std::cout << "garble: " << garble_out;
      //   std::cout << ", actual: " << hash_out << std::endl;
      // }

      assert(garble_out == hash_out);
    }
  }

  std::cout << "test hash passed" << std::endl;
}

void test_sort_complex(const int party) {
  const unsigned int mod_val = 513;
  Integer mod(32, mod_val, PUBLIC);

  const size_t size = 10;
  Integer *val_A = new Integer[size];
  Integer *val_B = new Integer[size];
  Integer *val = new Integer[size];
  Integer *freq_A = new Integer[size];
  Integer *freq_B = new Integer[size];
  Integer *freq = new Integer[size];

  // Input
  for (unsigned int i = 0; i < size; i++) {
    val_A[i] = Integer(32, rand()%mod_val, ALICE);
    val_B[i] = Integer(32, rand()%mod_val, BOB);
    freq_A[i] = Integer(32, rand()%100, ALICE);
    freq_B[i] = Integer(32, rand()%100, BOB);
  }

  // Compute
  for (unsigned int i = 0; i < size; i++) {
    val[i] = (val_A[i] + val_B[i]) % mod;
    freq[i] = (freq_A[i] + freq_B[i]) % mod;
  }
  sort(freq, size, val);

  // Reveal
  for(unsigned int i = 0; i < size; ++i) {
    int64_t v = val[i].reveal<int64_t>();
    int64_t f = freq[i].reveal<int64_t>();
    if (party == ALICE) {
      std::cout << "val[" << i << "] = " << v << ",\t freq = " << f << std::endl;
    }
  }
  delete[] val_A;
  delete[] val_B;
  delete[] val;
  delete[] freq_A;
  delete[] freq_B;
  delete[] freq;
}

void test_index(const int party) {
  const unsigned int mod_val = 513;
  const size_t nbits = LOG2(mod_val) + 1;  // 11
  Integer mod(nbits, mod_val, PUBLIC);

  const size_t size = 10;
  Integer* arr = new Integer[size];

  const unsigned int idx = 3;
  unsigned int actual = 0;

  // Input
  Integer index = Integer(nbits, idx, PUBLIC);

  for (unsigned int i = 0; i < size; i++) {
    int x = rand()%mod_val;
    arr[i] = Integer(nbits, x, PUBLIC);
    if (i == idx) actual = x;
  }

  // Get_at_index logic
  // Same mult gates as if have separate zero added instead.
  // I.e. res = If(b, arr, 0)
  Integer res(nbits, 0, PUBLIC);
  for (unsigned int i = 0; i < size; i++) {
    // nbits per
    Bit b = (index == Integer(nbits, i, PUBLIC));
    // nbits per
    res = If(b, arr[i], res);
  }

  // Reveal
  int64_t r = res.reveal<int64_t>();
  assert(actual == r);
  delete[] arr;

  std::cout << "Base index test passed" << std::endl;
}

void test_misc(const int party, int num_bits = 3) {

  // Turns out ints are kind of signed.

  int a_val = 8;
  int b_val = 7;

  Integer a(num_bits, a_val, PUBLIC);
  Integer b(num_bits, b_val, PUBLIC);
  Integer c = a % b;
  std::cout << num_bits << "-bit " << a_val << " % " << b_val << " = ";
  std::cout << c.reveal<int64_t>() << std::endl;
}

void test_full_sort(const int party, flint_rand_t hash_seed) {
  // Phase 0: Setup
  CountMinConfig cfg(.1, .4);  // 3 hashes range 7
  const size_t input_bits = 4;  // At least >= range bits
  CountMin count_min(cfg);
  count_min.init();
  count_min.setStore(input_bits, hash_seed);
  const size_t K = 3;

  // Phase 0.1: Populate count-min
  const size_t num_diff_vals = 5;
  unsigned int counts[num_diff_vals] = {1, 10, 5, 5, 2};
  int total_count = 0;
  for (unsigned int i = 0; i < num_diff_vals; i++) {
    count_min.add(i, counts[i]);
    total_count += counts[i];
  }

  // Shares in clear. Less garbled mod. About half as many gates
  const bool clear_version = false;

  /*
  Phase 0.2: Share split count-min and candidates
  Randomness is synced, so shares can be made locally.
  */
  fmpz_t share; fmpz_init(share);
  const unsigned int num_buckets = count_min.cfg.w * count_min.cfg.d;
  for (unsigned int i = 0; i < num_buckets; i++) {
    fmpz_randm(share, seed, Int_Modulus);
    if (party == ALICE) {
      fmpz_set(count_min.counts[i], share);
    } else {
      fmpz_mod_sub(count_min.counts[i], count_min.counts[i],
          share, mod_ctx);
    }
  }

  const size_t num_candidates = 10;
  int candidate_list[num_candidates] = {0, 0, 1, 2, 3, 1, 2, 3, 4, 4};
  // std::cout << "Creating " << num_candidates << " candidates " << std::endl;
  fmpz_t* candidates; new_fmpz_array(&candidates, num_candidates);
  for (unsigned int i = 0; i < num_candidates; i++) {
    // Set random. Turn off for clear
    if (!clear_version) fmpz_randm(candidates[i], seed, Int_Modulus);
    if (party == ALICE) {
      fmpz_set_ui(share, candidate_list[i]);
      fmpz_mod_sub(candidates[i], share, candidates[i], mod_ctx);
    }
  }
  fmpz_clear(share);

  // Run

  // Create struct
  HeavyEval heavy_eval(party, count_min, total_count);
  // heavy_eval.print_params(party == ALICE);

  // Step 1: de-share count-min
  heavy_eval.parse_countmin();
  // heavy_eval.print_countmin(party == ALICE);

  // Step 2: add values
  std::cout << "setting values..." << std::endl;
  if (clear_version) {
    heavy_eval.set_values_clear(candidates, num_candidates, ALICE);
  } else {
    heavy_eval.set_values(candidates, num_candidates);
  }
  // if (party == ALICE) std::cout << "parsed values: \n";
  // heavy_eval.print_values(party == ALICE);

  std::cout << "set values" << std::endl;

  // Step 3: Eval hashes
  const size_t method = 1;
  if (method == 0) {
    // Garble eval hashes
    heavy_eval.get_hashes();
  } else if (method == 1) {
    // Given secret shared candidates, eval hashes
    // Approximately 25% less garble
    //  - Method 0: Set value, eval hash, query freq
    //  - Method 1: Set value, set hash, query freq

    const size_t num_hashes = count_min.cfg.d;
    fmpz_t* hashes; new_fmpz_array(&hashes, num_candidates * num_hashes);
    fmpz_t x; fmpz_init(x);

    // Eval: count_min.store->eval(...);
    for (unsigned int i = 0; i < num_candidates; i++) {
      for (unsigned int j = 0; j < num_hashes; j++) {
        const size_t idx = i * num_hashes + j;
        count_min.store->eval(j, candidate_list[i], x);
        // Set random.
        if (!clear_version) fmpz_randm(hashes[idx], seed, Int_Modulus);
        // actual
        if (party == ALICE) {
          fmpz_mod_sub(hashes[idx], x, hashes[idx], mod_ctx);
        }
      }
    }

    if (clear_version) {
      heavy_eval.set_hashes_clear(hashes, ALICE);
    } else {
      heavy_eval.set_hashes(hashes);
    }

    clear_fmpz_array(hashes, num_candidates * num_hashes);
    fmpz_clear(x);
  }

  // Step 4: get frequencies
  heavy_eval.get_frequencies();
  // if (party == ALICE) std::cout << "parsed values with freqs:\n";
  // heavy_eval.print_values(party == ALICE);

  // Step 5: Sort and remove dupes
  heavy_eval.sort_remove_dupes(K);
  // if (party == ALICE) std::cout << "sorted values with freqs:\n";
  // heavy_eval.print_values(party == ALICE);

  uint64_t* top_values = new uint64_t[K];
  uint64_t* top_freqs = new uint64_t[K];

  // Step 6: extract out top K, return
  heavy_eval.return_top_K(K, top_values, top_freqs);

  // Check
  std::cout << "Top K = " << K << " values and freqs, decreasing\n";
  for (unsigned int i = 0; i < K; i++) {
    std::cout << "Value: " << top_values[i];
    std::cout << ", freq: " << top_freqs[i] << std::endl;
  }
  assert(top_values[0] == 1);
  assert(top_freqs[0] == 10);
  // 2 and 3 have same freq 5. Assert in some order
  assert(top_freqs[1] == 5);
  assert(top_freqs[2] == 5);
  assert((top_values[1] == 2 and top_values[2] == 3)
      or (top_values[1] == 3 and top_values[2] == 2));

  delete[] top_values;
  delete[] top_freqs;
  clear_fmpz_array(candidates, num_candidates);
}

void test_bucket_compare(const int party, flint_rand_t hash_seed) {
  const size_t input_bits = 4;
  const size_t output_range = 2;  // must be 2
  const size_t depth = input_bits;

  const size_t num_copies = 1;  // R
  const size_t num_groups = 4;  // < 2^input bits
  const size_t num_substreams = 2;
  const size_t num_pairs = num_copies * num_groups * num_substreams * depth;

  HashStoreBit store(num_groups, depth, input_bits, output_range, hash_seed);

  fmpz_t* bucket0; new_fmpz_array(&bucket0, num_pairs);
  fmpz_t* bucket1; new_fmpz_array(&bucket1, num_pairs);

  fmpz_t tmp; fmpz_init(tmp);
  fmpz_t tmp2; fmpz_init(tmp2);
  for (unsigned int i = 0; i < num_pairs; i++) {
    // actual[i] = abs(b0) < abs(b1);
    fmpz_randm(tmp, seed, Int_Modulus);
    fmpz_randm(tmp2, seed, Int_Modulus);
    if (party == ALICE) {
      fmpz_set(bucket0[i], tmp);
      fmpz_set(bucket1[i], tmp2);
    } else {
      // Actual
      int b0 = (i%2 ? (i+1) : (2*i+2)) * ((i/2)%2 ? 1 : -1);
      int b1 = (i%2 ? (2*i+2) : (i+1)) * (((i+1)/2)%2 ? 1 : -1);
      fmpz_set_si(bucket0[i], b0);
      fmpz_set_si(bucket1[i], b1);
      fmpz_mod(bucket0[i], bucket0[i], Int_Modulus);
      fmpz_mod(bucket1[i], bucket1[i], Int_Modulus);
      // Share adjust
      fmpz_mod_sub(bucket0[i], bucket0[i], tmp, mod_ctx);
      fmpz_mod_sub(bucket1[i], bucket1[i], tmp2, mod_ctx);
    }
  }

  HeavyExtract eval(party, store, num_copies, num_groups, num_substreams,
      depth, 2 * num_pairs + 2);
  // eval.print_params(party == ALICE);

  // Gates: 2 * (4 * share_bits - 2) * num_pairs
  eval.set_buckets(bucket0, bucket1);
  clear_fmpz_array(bucket0, num_pairs);
  clear_fmpz_array(bucket1, num_pairs);
  // eval.print_buckets(party == ALICE);

  // Gates: (6 * share_bits - 2 + freq_bits) * num_pairs
  eval.bucket_compare();
  // eval.print_cmp();

  // TODO: asserts

  fmpz_clear(tmp);
  fmpz_clear(tmp2);
}

void test_extract(const int party, const int sockfd, flint_rand_t hash_seed) {
  const size_t input_bits = 3;
  const size_t output_range = 2;  // test built only to support 2
  const size_t depth = input_bits;

  const size_t num_copies = 2;  // R
  const size_t num_groups = 4;  // < 2^input bits. Q
  const size_t num_substreams = 3;  // B
  const size_t num_candidates = num_copies * num_groups * num_substreams;
  const size_t num_pairs = num_candidates * depth;

  // Setup
  HashStoreBit store(num_groups, depth, input_bits, output_range, hash_seed);

  // figure out bits
  // For now, just replicate across copies, bits too small for split magic
  fmpz_t* bits; new_fmpz_array(&bits, num_pairs);
  fmpz_t tmp; fmpz_init(tmp);
  for (unsigned int i = 0; i < num_groups; i++) {
    for (unsigned int j = 0; j < depth; j++) {
      store.eval(i * depth + j, i, tmp);
      for (unsigned int k = 0; k < num_substreams; k++) {
        // Order group, substream, depth
        for (unsigned int r = 0; r < num_copies; r++) {
          const unsigned int idx = ((r * num_groups + i)
              * num_substreams + k) * depth + j;
          fmpz_set(bits[idx], tmp);
        }
      }
    }
  }

  // Set up bucket bits according to bits
  fmpz_t* bucket0; new_fmpz_array(&bucket0, num_pairs);
  fmpz_t* bucket1; new_fmpz_array(&bucket1, num_pairs);

  fmpz_t tmp2; fmpz_init(tmp2);
  for (unsigned int i = 0; i < num_pairs; i++) {
    fmpz_randm(tmp, seed, Int_Modulus);
    fmpz_randm(tmp2, seed, Int_Modulus);
    if (party == ALICE) {
      fmpz_set(bucket0[i], tmp);
      fmpz_set(bucket1[i], tmp2);
    } else {
      // Actuals
      // Base "random" value
      int b0 = (i/depth+1) * ((i/2)%2 ? 1 : -1);
      int b1 = (i/depth+1) * ((i/2)%2 ? 1 : -1);
      // Set actual based on desired bits
      (fmpz_is_one(bits[i]) ? b1 : b0) *= 2;
      fmpz_set_si(bucket0[i], b0);
      fmpz_set_si(bucket1[i], b1);
      fmpz_mod(bucket0[i], bucket0[i], Int_Modulus);
      fmpz_mod(bucket1[i], bucket1[i], Int_Modulus);
      // Share adjust
      fmpz_mod_sub(bucket0[i], bucket0[i], tmp, mod_ctx);
      fmpz_mod_sub(bucket1[i], bucket1[i], tmp2, mod_ctx);
    }
  }
  clear_fmpz_array(bits, num_pairs);
  fmpz_clear(tmp);
  fmpz_clear(tmp2);

  // 2 on one side, 1 on the other, for +3
  // n(n+1)/2, with i+1 stuff
  const unsigned int total = 3 * (num_candidates + 1) * num_candidates / 2;

  const int method = 1;

  if (method == 0) {
    // Run
    HeavyExtract eval(party, store, num_copies, num_groups, num_substreams,
        depth, total);
    // eval.print_params(party == ALICE);

    eval.set_buckets(bucket0, bucket1);
    clear_fmpz_array(bucket0, num_pairs);
    clear_fmpz_array(bucket1, num_pairs);
    // eval.print_buckets(party == ALICE);

    eval.bucket_compare();
    // eval.print_cmp(party == ALICE);

    eval.extract_candidates();
    // eval.print_candidates(party == ALICE);

    Integer* candidates = eval.get_candidates();
    for (unsigned int i = 0; i < num_candidates; i++) {
      uint64_t v = candidates[i].reveal<uint64_t>();
      uint64_t expected = (i/num_substreams) % num_groups;
      // std::cout << "got[" << i << "]: " << v;
      // std::cout << ", expected = " << expected << std::endl;
      assert(v == expected);
    }
  } else if (method == 1) {
    bool* cmp = bucket_compare_clear(sockfd, num_pairs, bucket0, bucket1);
    clear_fmpz_array(bucket0, num_pairs);
    clear_fmpz_array(bucket1, num_pairs);
    uint64_t* candidates = new uint64_t[num_candidates];
    extract_candidates_clear(num_copies, num_groups, num_substreams, depth,
        input_bits, store, cmp, candidates);
    delete[] cmp;

    for (unsigned int i = 0; i < num_candidates; i++) {
      uint64_t v = candidates[i];
      uint64_t expected = (i/num_substreams) % num_groups;
      // std::cout << "got[" << i << "]: " << v;
      // std:: << ", expected = " << expected << std::endl;
      assert(v == expected);
    }

    delete[] candidates;
  }
}

void run(const int party) {
  init_constants();
  // Fork. Party 1, 2. same port

  flint_rand_t hash_seed;
  flint_randinit(hash_seed);

  int sockfd = 0; int tmp_sock = 0;
  if (party == ALICE) {
    std::cout << "ALICE" << std::endl;
    tmp_sock = init_receiver();
    sockfd = accept_receiver(tmp_sock);
  } else {
    std::cout << "BOB" << std::endl;
    sockfd = init_sender(SERVER0_IP);
  }

  OT_Wrapper* ot0 = new OT_Wrapper(party == ALICE ? nullptr : SERVER0_IP, 60051);
  OT_Wrapper* ot1 = new OT_Wrapper(party == BOB ? nullptr : SERVER1_IP, 60052);

  // Somehow new io, and ot1->io dont' like to work. Grab from OT
  NetIO* io = ot0->io;

  setup_semi_honest(io, party);

  // Use some OT, since garble shares IO with OT.
  // So we ensure it still works being shared
  uint64_t data0[2] = {10, 20};
  uint64_t data0_1[2] = {11, 21};
  uint64_t data1[2] = {1000, 2000};
  uint64_t data1_1[2] = {1001, 2001};
  uint64_t* data = new uint64_t[2];
  uint64_t* data_1 = new uint64_t[2];
  bool c[2] = {false, true};
  if (party == ALICE) {
    ot0->send(data0, data1, 2, data0_1, data1_1);
    ot1->recv(data, c, 2, data_1);
  } else {
    ot0->recv(data, c, 2, data_1);
    ot1->send(data0, data1, 2, data0_1, data1_1);
  }
  delete[] data;
  delete[] data_1;

  /* General Behavior */
  // test_compare(party);
  // test_mod(party);
  // test_sort_complex(party);
  // test_index(party);
  // test_misc(party, 3);
  // test_misc(party, 4);
  // test_misc(party, 5);
  // test_misc(party, 6);

  /* EvalHeavy specific */
  // test_hash(party, hash_seed);
  test_full_sort(party, hash_seed);

  io->flush();  // Need to flush before using sockets
  // test_bucket_compare(party, hash_seed);
  test_extract(party, sockfd, hash_seed);

  io->flush();

  std::cout << "Mult gates: " << CircuitExecution::circ_exec->num_and() << std::endl;

  finalize_semi_honest();  // just deletes things
  // delete io;  // gets deleted with ot0
  delete ot0;
  delete ot1;
  close(sockfd); if (party == ALICE) close(tmp_sock);
  clear_constants();
}


int main(int argc, char** argv) {
  if (argc == 1) {
    // Note: Garble doesn't work well with threading, it seems
    pid_t pid = fork();
    if (pid == 0) {
      run(ALICE);
    } else {
      run(BOB);
    }

    return 0;
  }

  int party = atoi(argv[1]) + 1;
  std::cout << "Party: " << party << std::endl;
  run(party);
}
