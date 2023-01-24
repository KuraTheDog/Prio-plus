#undef NDEBUG
#include <assert.h>
#include <iostream>

#include "../fmpz_utils.h"
#include "../hash.h"

const size_t input_bits = 3;
const size_t output_range = 5;
const size_t num_hashes = 3;

void run_hash(HashStore &store, unsigned int i) {
  fmpz_t out; fmpz_init(out);
  unsigned int max = (1ULL << input_bits);
  if (max > 8) max = 8;
  for (unsigned int j = 0; j < max; j++) {
    store.eval(i, j, out);
    std::cout << "Eval_" << i << "(" << j << ") \t= ";
    fmpz_print(out);
    std::cout << std::endl;
  }
  fmpz_clear(out);
}

void test_seed_sync() {
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

  std::cout << "First store" << std::endl;
  for (unsigned int i = 0; i < num_hashes; i++) {
    store.print_hash(i);
    run_hash(store, i);
  }
  // store.print();

  // std::cout << "Second store" << std::endl;
  // store2.print();
}

// void test_inverse() {
//   flint_rand_t hash_seed;
//   flint_randinit(hash_seed);

//   HashStoreBit store(d, l_bits, w, hash_seed);
// }

int main(int argc, char** argv){
  test_seed_sync();

  // test_inverse();
}
