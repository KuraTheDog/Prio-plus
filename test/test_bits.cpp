#include <iostream>

#include "../constants.h"
#include "../edabit.h"
#include "../server.h"

void localTest() {
int n = 6;
  int pow = 1 << n;
  std::cout << "n = " << n << ", so pow = " << pow << std::endl;

  fmpz_t tmp;
  fmpz_init(tmp);

  // random adjusting. comment to change seed.
  // fmpz_randm(tmp, seed, Int_Modulus);
  // fmpz_randm(tmp, seed, Int_Modulus);

  EdaBit* ebit00 = new EdaBit(n);  // 0's share of 0
  EdaBit* ebit01 = new EdaBit(n);  // 1's share of 0
  makeLocalEdaBit(ebit00, ebit01, n);

  EdaBit* ebit10 = new EdaBit(n);  // 0's share of 1
  EdaBit* ebit11 = new EdaBit(n);  // 1's share of 1
  makeLocalEdaBit(ebit10, ebit11, n);

  // 0 has ebit00, ebit01. 1 has ebit10, ebit11
  // They move things around
  // 0 has ebit00, ebit10. 1 has ebit01, ebit11

  std::cout << "0 has ebit00: \n"; ebit00->print();
  std::cout << "0 has ebit10: \n"; ebit10->print();

  std::cout << "1 has ebit01: \n"; ebit01->print();
  std::cout << "1 has ebit11: \n"; ebit11->print();

  fmpz_add(tmp, ebit00->r, ebit01->r);
  fmpz_mod(tmp, tmp, Int_Modulus);
  std::cout << "ebit00 and ebit01, sum = "; fmpz_print(tmp); std::cout << std::endl;
  std::cout << "ebit00 and ebit01, xor = " << (ebit00->get_int_b() ^ ebit01->get_int_b()) << std::endl;

  fmpz_add(tmp, ebit10->r, ebit11->r);
  fmpz_mod(tmp, tmp, Int_Modulus);
  std::cout << "ebit10 and ebit11, sum = "; fmpz_print(tmp); std::cout << std::endl;
  std::cout << "ebit10 and ebit11, xor = " << (ebit10->get_int_b() ^ ebit11->get_int_b()) << std::endl;

  // Final bits
  EdaBit* ebit0 = new EdaBit(n);
  EdaBit* ebit1 = new EdaBit(n);

  // r'
  fmpz_add(ebit0->r, ebit00->r, ebit10->r);
  fmpz_mod(ebit0->r, ebit0->r, Int_Modulus);
  std::cout << "[r']^0 = "; fmpz_print(ebit0->r); std::cout << std::endl;

  fmpz_add(ebit1->r, ebit01->r, ebit11->r);
  fmpz_mod(ebit1->r, ebit1->r, Int_Modulus);
  std::cout << "[r']^1 = "; fmpz_print(ebit1->r); std::cout << std::endl;

  // bs
  // x0 = ebit00->b, y0 = ebit10->b
  // x1 = ebit00->b, y1 = ebit10->b
  // z0 = ebit0->b, z1 = ebit1->b
  // c_{i+1} = c_i xor ((x_i xor c_i) and (y_i xor c_i))
  // output z_i = x_i xor y_i xor c_i

  bool carry0 = false, carry1 = false;  // carry bits
  bool x0, x1, y0, y1, tmp0, tmp1;  // addition calculation
  bool d, e, d0, d1, e0, e1;  // triple calculation
  bool a0, a1, b0, b1, c0, c1;  // beaver triple;

  for (int i = 0; i < n; i++) {
    ebit0->b[i] = carry0 ^ ebit00->b[i] ^ ebit10->b[i];
    x0 = carry0 ^ ebit00->b[i];
    y0 = carry0 ^ ebit10->b[i];

    ebit1->b[i] = carry1 ^ ebit01->b[i] ^ ebit11->b[i];
    x1 = carry1 ^ ebit01->b[i];
    y1 = carry1 ^ ebit11->b[i];

    // boolean beaver triple.
    // (a0^a1) * (b0^b1) = (c0^c1). 5 random per.
    fmpz_randbits(tmp, seed, 5);
    a0 = get_fmpz_bit(tmp, 0);
    a1 = get_fmpz_bit(tmp, 1);
    b0 = get_fmpz_bit(tmp, 2);
    b1 = get_fmpz_bit(tmp, 3);
    c0 = get_fmpz_bit(tmp, 4);
    c1 = c0 ^ ((a0 ^ a1) and (b0 ^ b1));

    // [c] xor [x]*e xor [y] * d xor e * d
    
    d0 = x0 ^ a0;
    e0 = y0 ^ b0;

    d1 = x1 ^ a1;
    e1 = y1 ^ b1;
    // broadcast
    d = d0 ^ d1;
    e = e0 ^ e1;

    tmp0 = c0 ^ (x0 and e) ^ (y0 and d) ^ (d and e);
    carry0 ^= tmp0;

    tmp1 = c1 ^ (x1 and e) ^ (y1 and d);
    carry1 ^= tmp1;
  }

  // TODO: Magical adder circuit, rather than add.
  std::cout << " [b]_2^0 = " << ebit0->get_int_b() << (carry0? " with carry" : " no carry") << std::endl;
  std::cout << " [b]_2^1 = " << ebit1->get_int_b() << (carry1? " with carry" : " no carry") << std::endl;

  // b_p's

  // dabit for adjusting. Doing one bit convertToField
  // TODO: actually do via servers
  DaBit* dabit0 = new DaBit();
  DaBit* dabit1 = new DaBit();
  makeLocalDaBit(dabit0, dabit1);
  std::cout << "conversion dabit0\n"; dabit0->print();
  std::cout << "conversion dabit1\n"; dabit1->print();
  bool v0 = carry0 ^ dabit0->b2;
  bool v1 = carry1 ^ dabit1->b2;
  // 0 has v0, 1 has v1. Broadcast both so both know v.
  bool v = v0 ^ v1;
  fmpz_t bp0;
  if (v) {  // v = 1, so x0 = 1 - [b]_p
    fmpz_init_set_si(bp0, 1);
    fmpz_sub(bp0, bp0, dabit0->bp);
    fmpz_mod(bp0, bp0, Int_Modulus);
  } else {  // v = 0, so x0 = [b]_p
    fmpz_init_set(bp0, dabit0->bp);
  }
  std::cout << "bp0 = "; fmpz_print(bp0); std::cout << std::endl;
  fmpz_t bp1;
  if (v) {  // v = 1, so x1 = [b]_p
    fmpz_init_set_si(bp1, 0);
    fmpz_sub(bp1, bp1, dabit1->bp);
    fmpz_mod(bp1, bp1, Int_Modulus);
  } else {  // v = 0, so x1 = [b]_p
    fmpz_init_set(bp1, dabit1->bp);
  }
  std::cout << "bp1 = "; fmpz_print(bp1); std::cout << std::endl;

  // r = r' - 2^n * carry * bp
  if (carry0) {
    fmpz_mul_ui(tmp, bp0, pow);
    fmpz_sub(ebit0->r, ebit0->r, tmp);
    fmpz_mod(ebit0->r, ebit0->r, Int_Modulus);
  }

  if (carry1) {
    fmpz_mul_ui(tmp, bp1, pow);
    fmpz_sub(ebit1->r, ebit1->r, tmp);
    fmpz_mod(ebit1->r, ebit1->r, Int_Modulus);
  }

  std::cout << "Final ebit0:\n"; ebit0->print();
  std::cout << "Final ebit1:\n"; ebit1->print();

  fmpz_add(tmp, ebit0->r, ebit1->r);
  fmpz_mod(tmp, tmp, Int_Modulus);
  std::cout << "ebit0 and ebit1, sum = "; fmpz_print(tmp); std::cout << std::endl;
  std::cout << "ebit0 and ebit1, xor = " << (ebit0->get_int_b() ^ ebit1->get_int_b()) << std::endl;
}

int main(int argc, char* argv[]) {
  init_constants();
  
  localTest();

  clear_constants();
  return 0;
}