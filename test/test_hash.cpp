#undef NDEBUG
#include <assert.h>
#include <iostream>

#include "../fmpz_utils.h"
#include "../hash.h"
#include "../heavy.h"
#include "../zipf.h"

const size_t input_bits = 3;
const size_t output_range = 5;

void run_hash(HashStore &store, unsigned int i) {
  fmpz_t out; fmpz_init(out);
  unsigned int max = (1ULL << input_bits);
  if (max > 8) max = 8;
  for (unsigned int j = 0; j < max; j++) {
    store.eval(i, j, out);
    std::cout << "Eval_" << i << "(" << j << ") \t= " << fmpz_get_ui(out) << "\n";
  }
  fmpz_clear(out);
}

void test_seed_sync() {
  const size_t num_hashes = 12;

  flint_rand_t hash_seed;
  flint_randinit(hash_seed);

  // rand adjust
  fmpz_t tmp; fmpz_init(tmp);
  for (unsigned int i = 0; i < 5; i++)
    fmpz_randbits(tmp, hash_seed, 1);

  // Ensure seed clone works: aka both servers have seed
  flint_rand_t second_hash_seed;
  flint_randinit(second_hash_seed);
  second_hash_seed[0] = hash_seed[0];

  HashStorePoly store(num_hashes, input_bits, output_range, hash_seed);
  HashStorePoly store2(num_hashes, input_bits, output_range, second_hash_seed);

  fmpz_t out1; fmpz_init(out1);
  fmpz_t out2; fmpz_init(out2);

  store.eval(1, 2, out1);
  store2.eval(1, 2, out2);
  assert(fmpz_equal(out1, out2));

  // std::cout << "First store" << std::endl;
  // for (unsigned int i = 0; i < num_hashes; i++) {
  //   store.print_hash(i);
  //   run_hash(store, i);
  // }
  // store.print();

  // std::cout << "Second store" << std::endl;
  // store2.print();
}

void test_inverse(size_t group_size) {
  group_size = group_size ? group_size : input_bits;
  flint_rand_t hash_seed;
  flint_randinit(hash_seed);
  fmpz_t tmp; fmpz_init(tmp);

  // Offset randomness
  for (unsigned int i = 0; i < 5; i++)
    fmpz_randbits(tmp, hash_seed, 1);
  fmpz_clear(tmp);

  const size_t num_groups = 3;  // < 2^input bits, for testing

  std::cout << "testing inverse with group size = " << group_size << std::endl;
  HashStoreBit store(num_groups, group_size, input_bits, output_range, hash_seed);

  // std::cout << "Coeff matrix: " << std::endl;
  // store.print_coeff();

  fmpz_t* values; new_fmpz_array(&values, group_size);
  for (unsigned int i = 0; i < num_groups; i++) {
    // std::cout << "Testing group: " << i << std::endl;
    const unsigned int x = i + 1;
    // std::cout << "input is " << x << std::endl;
    for (unsigned int j = 0; j < group_size; j++) {
      store.eval(i * group_size + j, x, values[j]);
      // store.print_hash(j);
      // std::cout << "eval_" << j << "(" << x << ") = ";
      // std::cout << fmpz_get_ui(values[j]) << std::endl;
    }

    uint64_t ans;
    int ret = store.solve(i, values, ans);
    // std::cout << "Answer is " << ans << ", with " << ret << " invalid" << std::endl;
    assert(ret == 0);
    assert(ans == x);
  }
  clear_fmpz_array(values, group_size);
}

void test_countmin() {
  flint_rand_t hash_seed;
  flint_randinit(hash_seed);

  // .1, .1 -> 3 hash range 28
  // So input bits becomes 5
  CountMinConfig cfg(0.1, 0.1);
  cfg.print();

  CountMin count(cfg);
  count.init();
  count.setStore(3, hash_seed);

  // count.store->print();

  // Small (1/w)^d chance of each full overlapping. so test not too large to amplify this.
  uint64_t vals[3] = {1, 2, 3};
  unsigned int counts[3] = {6, 2, 0};

  for (unsigned int i = 0; i < 3; i++)
    count.add(vals[i], counts[i]);

  count.print();

  for (unsigned int i = 0; i < 3; i++) {
    unsigned int ans = count.query(vals[i]);
    // std::cout << "query(" << vals[i] << ") = " << ans;
    // std::cout  << ", vs actual " << counts[i] << std::endl;
    assert(ans == counts[i]);
  }
}

void test_half() {
  flint_rand_t hash_seed;
  flint_randinit(hash_seed);

  const size_t num_bits = 3;

  HashStoreHalf store(num_bits, hash_seed);

  for (unsigned int i = 0; i < num_bits+1; i++) {
    unsigned int n = 0;
    for (unsigned int x = 0; x < (1 << num_bits); x++)
      n += (unsigned int) store.survives(i, x);
    // Right number of survivors. Odd ax+b forces to be true.
    // std::cout << i << ", " << n << std::endl;
    assert(n == 1ULL << (num_bits - i));
  }
}

void test_zipf() {
  const int support = 10000;
  const double exponent = 1.03;

  ZipF distribution(support, exponent);

  const int N = 1000000;
  const int K = 15;
  int count[K];
  memset(count, 0, K * sizeof(int));

  std::cout << "Zipf(" << support << ", " << exponent << ")\n";

  for (unsigned int i = 0; i < N; i++) {
    int val = distribution.sample(5 * K);
    if (val < K) count[val]++;
  }
  for (unsigned int i = 0; i < K; i++) {
    std::cout << "Zipf[" << i << "] = " << count[i];
    std::cout << " -> " << 1. * count[i] / N << "\n";
  }
}

int main(int argc, char** argv){
  test_seed_sync();

  test_countmin();
  test_inverse(6);
  test_inverse(0);

  test_half();

  test_zipf();
}
