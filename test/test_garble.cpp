#undef NDEBUG
#include <iostream>

#include "constants.h"
#include "eval_heavy.h"
#include "hash.h"
#include "heavy.h"

#include "emp-sh2pc/emp-sh2pc.h"

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

void test_compare(int party) {

  Integer a(4, 3, ALICE);
  Integer b(4, 2, BOB);
  Bit res = a > b;

  std::cout << "Party " << party << ", ALICE larger?\t"<< res.reveal<bool>()<<endl;
}

void test_hash(int party, flint_rand_t hash_seed) {

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
      //   std::cout << "garble: " << garble_out << ", actual: " << hash_out << std::endl;
      // }

      assert(garble_out == hash_out);
    }
  }

  std::cout << "test hash passed" << std::endl;
}

void test_sort_complex(int party) {
  int mod_val = 513;
  Integer mod(32, mod_val, PUBLIC);

  int size = 10;
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
  for(int i = 0; i < size; ++i) {
    int64_t v = val[i].reveal<int64_t>();
    int64_t f = freq[i].reveal<int64_t>();
    if (party == ALICE) {
      std::cout << "val[" << i << "] = " << v << ",\t freq = " << f << std::endl;
    }
  }
}

void test_index(int party) {
  int mod_val = 513;
  Integer mod(32, mod_val, PUBLIC);

  int size = 10;
  Integer* arr_A = new Integer[size];
  Integer* arr_B = new Integer[size];
  Integer* arr = new Integer[size];

  int idx = 3;
  int actual;

  // Input
  Integer index = Integer(32, idx, PUBLIC);

  for (unsigned int i = 0; i < size; i++) {
    int a = rand()%mod_val;
    int b = rand()%mod_val;
    arr_A[i] = Integer(32, a, ALICE);
    arr_B[i] = Integer(32, b, BOB);
    if (i == idx) actual = (a + b) % mod_val;
  }

  // Compute
  for (unsigned int i = 0; i < size; i++) {
    arr[i] = (arr_A[i] + arr_B[i]) % mod;
  }

  // Same mult gates as if have separate zero added instead.
  // I.e. res = If(b, arr, 0)
  Integer res(32, 0, PUBLIC);
  for (unsigned int i = 0; i < size; i++) {
    Bit b = (index == Integer(32, i, PUBLIC));

    res = If(b, arr[i], res);
  }

  // Reveal
  int64_t r = res.reveal<int64_t>();
  // if (party == ALICE) {
  //   std::cout << "Mult gates: " << CircuitExecution::circ_exec->num_and()<<endl;
  //   std::cout << "result at index " << idx << " = " << r << std::endl;
  // }
  // for(int i = 0; i < size; ++i) {
  //   int64_t a = arr[i].reveal<int64_t>();
  //   if (party == ALICE) {
  //     std::cout << "arr[" << i << "] = " << a << std::endl;
  //   }
  // }
  assert(actual == r);
  delete[] arr_A;
  delete[] arr_B;
  delete[] arr;

  std::cout << "Base index test passed" << std::endl;
}

void test_misc(int party, int num_bits = 3) {

  // Turns out ints are kind of signed.

  int a_val = 8;
  int b_val = 7;

  Integer a(num_bits, a_val, PUBLIC);
  Integer b(num_bits, b_val, PUBLIC);
  Integer c = a % b;
  std::cout << num_bits << "-bit " << a_val << " % " << b_val << " = ";
  std::cout << c.reveal<int64_t>() << std::endl;
}

void test_full(int party, flint_rand_t hash_seed) {

  // Phase 0: Setup
  CountMinConfig cfg(.1, .4);  // 3 hashes range 7
  const size_t input_bits = 4;  // At least >= range bits
  CountMin count_min(cfg);
  count_min.init();
  count_min.setStore(input_bits, hash_seed);

  // Phase 0.1: Populate count-min
  int vals[3] = {1, 2, 3};
  unsigned int counts[3] = {1, 6, 3};
  int total_count = 0;
  for (unsigned int i = 0; i < 3; i++) {
    count_min.add(vals[i], counts[i]);
    total_count += counts[i];
  }

  /*
  Phase 0.2: Share split count-min and candidates
  Randomness is synced, so shares can be made locally.
  */
  fmpz_t share; fmpz_init(share);
  const int num_buckets = count_min.cfg.w * count_min.cfg.d;
  for (unsigned int i = 0; i < num_buckets; i++) {
    fmpz_randm(share, seed, Int_Modulus);
    if (party == ALICE) {
      fmpz_set(count_min.counts[i], share);
    } else {
      fmpz_mod_sub(count_min.counts[i], count_min.counts[i],
          share, mod_ctx);
    }
  }
  // Candidates are 0, 1, 2, 3, ...
  const size_t num_candidates = 8;
  std::cout << "Creating " << num_candidates << "c andidates " << std::endl;
  fmpz_t* candidates; new_fmpz_array(&candidates, num_candidates);
  for (unsigned int i = 0; i < num_candidates; i++) {
    fmpz_randm(share, seed, Int_Modulus);
    if (party == ALICE) {
      fmpz_set(candidates[i], share);
    } else {
      int v = (i/2) + 1;
      fmpz_set_ui(candidates[i], v);
      fmpz_mod_sub(candidates[i], candidates[i], share, mod_ctx);
    }
  }

  // Create struct
  HeavyEval heavy_eval(party, count_min, total_count);

  // Step 1: de-share count-min
  heavy_eval.parse_countmin();
  // heavy_eval.print_countmin();

  // Step 2: add values
  heavy_eval.set_values(candidates, num_candidates);
  // std::cout << "parsed values: \n";
  // heavy_eval.print_values();

  // Step 3: get frequencies
  heavy_eval.get_frequencies();
  // std::cout << "parsed values with freqs: \n";
  // heavy_eval.print_values();

  heavy_eval.sort_remove_dupes();
  // std::cout << "sorted values with freqs: \n";
  // heavy_eval.print_values();

  const size_t K = 3;
  uint64_t* top_values = new uint64_t[K];
  uint64_t* top_freqs = new uint64_t[K];
  heavy_eval.return_top_K(K, top_values, top_freqs);
  if (party == ALICE) {
    std::cout << "Top K = " << K << " values and freqs, decreasing\n";
    for (unsigned int i = 0; i < K; i++) {
      std::cout << "Value: " << top_values[i] << ", freq: " << top_freqs[i] << std::endl;
    }
  }

  clear_fmpz_array(candidates, num_candidates);
  fmpz_clear(share);
}


int main(int argc, char** argv) {
  init_constants();
  int port = 12345;
  int party;
  // Fork. Party 1, 2. same port. Not using sockets

  flint_rand_t hash_seed;
  flint_randinit(hash_seed);

  pid_t pid = fork();

  if (pid == 0) {
    party = emp::ALICE;
  } else {
    party = emp::BOB;
  }

  NetIO* io = new NetIO(party==ALICE ? nullptr : "127.0.0.1", port);
  setup_semi_honest(io, party);

  /* General Behavior */
  // test_compare(party);
  // test_sort_complex(party);
  // test_index(party);
  // test_misc(party, 3);
  // test_misc(party, 4);
  // test_misc(party, 5);
  // test_misc(party, 6);

  /* EvalHeavy specific */
  // test_hash(party, hash_seed);
  test_full(party, hash_seed);

  if (party == ALICE) {
    std::cout << "Mult gates: " << CircuitExecution::circ_exec->num_and()<<endl;
  }

  finalize_semi_honest();
  delete io;
  clear_constants();
}
