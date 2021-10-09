#include <iostream>

#include "../fmpz_utils.h"
#include "../hash.h"

const size_t l_bits = 4;
const size_t w = 5;
const size_t d = 3;

void run_hash(HashStore &store, unsigned int i) {
  fmpz_t out; fmpz_init(out);
  unsigned int max = (1ULL << l_bits);
  if (max > 8) max = 8;
  for (unsigned int j = 0; j < max; j++) {
    store.eval(i, j, out);
    std::cout << "Eval(" << i << ", " << j << ") \t= ";
    fmpz_print(out);
    std::cout << std::endl;
  }
  fmpz_clear(out);
}

void test_hash() {
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

  std::cout << "First store" << std::endl;
  HashStore store(d, l_bits, w, hash_seed);
  for (unsigned int i = 0; i < d; i++) {
    store.print_hash(i);
    run_hash(store, i);
  }

  std::cout << "Second store" << std::endl;
  HashStore store2(d, l_bits, w, second_hash_seed);
  for (unsigned int i = 0; i < d; i++)
    store.print_hash(i);
}

int main(int argc, char** argv){
  test_hash();
}
