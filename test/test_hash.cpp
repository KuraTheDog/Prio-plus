#include <iostream>

#include "../fmpz_utils.h"
#include "../hash.h"

const size_t w = 8;
const size_t d = 4;

void run_hash(HashStore &store, unsigned int i) {
  fmpz_t x; fmpz_init(x);
  fmpz_t out; fmpz_init(out);
  for (unsigned int j = 0; j < 3; j++) {
    fmpz_set_ui(x, j);
    store.eval(i, x, out);
    std::cout << "Eval(" << i << ", " << j << ") = ";
    fmpz_print(out);
    std::cout << std::endl;
  }
  fmpz_clear(x);
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
  HashStore store(w, d, hash_seed);
  for (unsigned int i = 0; i < d; i++)
    store.print_hash(i);

  std::cout << "Second store" << std::endl;
  HashStore store2(w, d, second_hash_seed);
  for (unsigned int i = 0; i < d; i++)
    store.print_hash(i);

  flint_randclear(hash_seed);
  flint_randclear(second_hash_seed);
}

int main(int argc, char** argv){
  test_hash();
}