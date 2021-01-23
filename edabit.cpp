/* edaBit related logic, for share conversion.
Based on ia.cr/2020/338
*/

#include "constants.h"
#include "edabit.h"
#include "net_share.h"

bool multiplyBoolShares(const int serverfd, const int server_num, const bool x, const bool y, const BooleanBeaverTriple* triple) {
  bool z, d, e, d_this, e_this, d_other, e_other;

  d_this = x ^ triple->a;
  e_this = y ^ triple->b;

  send_bool(serverfd, d_this);
  recv_bool(serverfd, d_other);
  send_bool(serverfd, e_this);
  recv_bool(serverfd, e_other);

  d = d_this ^ d_other;
  e = e_this ^ e_other;

  z = triple->c ^ (x and e) ^ (y and d);
  if (server_num == 0)
      z ^= (d and e);

  return z;
}

void multiplyArithmeticShares(const int serverfd, const int server_num, const fmpz_t x, const fmpz_t y, fmpz_t z, const BeaverTriple* triple) {
  fmpz_t d, d_other, e, e_other;

  fmpz_init(d);
  fmpz_init(e);
  fmpz_init(d_other);
  fmpz_init(e_other);

  fmpz_add(d, x, triple->A);  // [x] - [a]
  fmpz_mod(d, d, Int_Modulus);
  fmpz_add(e, y, triple->B);  // [y] - [b]
  fmpz_mod(e, e, Int_Modulus);

  send_fmpz(serverfd, d);
  recv_fmpz(serverfd, d_other);
  send_fmpz(serverfd, e);
  recv_fmpz(serverfd, e_other);

  fmpz_add(d, d, d_other);  // x - a
  fmpz_mod(d, d, Int_Modulus);
  fmpz_add(e, e, e_other);  // y - b
  fmpz_mod(e, e, Int_Modulus);

  // [xy] = [c] + [x] e + [y] d - de
  fmpz_set(z, triple->C);
  fmpz_addmul(z, x, e);
  fmpz_addmul(z, y, d);
  if (server_num == 0)
      fmpz_submul(z, d, e);
  fmpz_mod(z, z, Int_Modulus);

  fmpz_clear(d);
  fmpz_clear(d_other);
  fmpz_clear(e);
  fmpz_clear(e_other);
}

// c_{i+1} = c_i xor ((x_i xor c_i) and (y_i xor c_i))
// output z_i = x_i xor y_i xor c_i
bool addBinaryShares(const int serverfd, const int server_num, const size_t n, const bool* x, const bool* y, bool* z, const BooleanBeaverTriple* triples) {
  bool carry = false;
  for (int i = 0; i < n; i++) {
    z[i] = carry ^ x[i] ^ y[i];
    bool xi = carry ^ x[i];
    bool yi = carry ^ y[i];
    carry ^= multiplyBoolShares(serverfd, server_num, xi, yi, &triples[i]);
  }
  return carry;
}

void b2a_daBit(const int serverfd, const int server_num, const DaBit* dabit, const bool x, fmpz_t &xp) {
  const bool v_this = x ^ dabit->b2;
  bool v_other;
  send_bool(serverfd, v_this);
  recv_bool(serverfd, v_other);
  const bool v = v_this ^ v_other;

  // [x]_p = v + [b]_p - 2 v [b]_p. v only added for one server.
  // So we only add v on server 1 (when v = 1)
  if (v) {  // If v = 1, then [x]_p = (0/1) - [b]_p
    fmpz_set_ui(xp, server_num);
    fmpz_sub(xp, xp, dabit->bp);
    fmpz_mod(xp, xp, Int_Modulus);
  } else {  // If v = 0, then [x]_p = [b]_p
    fmpz_set(xp, dabit->bp);
  }
}

void b2a_edaBit(const int serverfd, const int server_num, const EdaBit* edabit, const fmpz_t x, fmpz_t &xp, const BooleanBeaverTriple* triples) {
  const size_t n = edabit->n;

  // Convert x2 to bool array
  bool x2[n];
  for (int i = 0; i < n; i++)
    x2[i] = get_fmpz_bit(x, i);

  // [x+r]_2 = [x]_2 + [r]_2 via circuit
  bool xr[n+1];
  bool carry = addBinaryShares(serverfd, server_num, n, x2, edabit->b, xr, triples);
  xr[n] = carry;

  fmpz_from_bool_array(xp, xr, n+1);  // [x + r]_2

  // Reveal x + r, convert to mod p shares
  if (server_num == 0) {
    fmpz_t xr_other;
    fmpz_init(xr_other);
    recv_fmpz(serverfd, xr_other);  // get other [x + r]_2
    fmpz_xor(xp, xp, xr_other);     // Real x + r
    fmpz_randm(xr_other, seed, Int_Modulus);  // other [x + r]_p
    fmpz_sub(xp, xp, xr_other);
    fmpz_mod(xp, xp, Int_Modulus);  // This [x + r]_p
    send_fmpz(serverfd, xr_other);
    fmpz_clear(xr_other);
  } else {
    send_fmpz(serverfd, xp);
    recv_fmpz(serverfd, xp);  // This [x + r]_p
  }

  // [x]_p = [x+r]_p - [r]_p
  fmpz_sub(xp, xp, edabit->r);
  fmpz_mod(xp, xp, Int_Modulus);
}

DaBit* generateDaBit(const int serverfd, const int server_num, const BeaverTriple* triple) {
  // Answer
  DaBit* dabit = new DaBit();

  // Create local
  DaBit* dabit0 = new DaBit();
  DaBit* dabit1 = new DaBit();
  makeLocalDaBit(dabit0, dabit1);

  // Exchange
  DaBit* tmp_dabit = new DaBit();
  if (server_num == 0) {
    send_DaBit(serverfd, dabit1);
    recv_DaBit(serverfd, tmp_dabit);
    delete dabit1;
    dabit1 = tmp_dabit;
  } else {
    send_DaBit(serverfd, dabit0);
    recv_DaBit(serverfd, tmp_dabit);
    delete dabit0;
    dabit0 = tmp_dabit;
  }

  // Xor boolean shares
  dabit->b2 = dabit0->b2 ^ dabit1->b2;

  // Xor arithmetic shares, using a xor b = a + b - 2ab
  fmpz_t z;
  fmpz_init(z);
  multiplyArithmeticShares(serverfd, server_num, dabit0->bp, dabit1->bp, z, triple);

  fmpz_add(dabit->bp, dabit0->bp, dabit1->bp);
  fmpz_submul_ui(dabit->bp, z, 2);
  fmpz_mod(dabit->bp, dabit->bp, Int_Modulus);

  fmpz_clear(z);
  delete dabit0;
  delete dabit1;

  return dabit;
}

EdaBit* generateEdaBit(const int serverfd, const int server_num, const size_t n, const BooleanBeaverTriple* triples, const DaBit* dabit) {
  // Answer bit
  EdaBit* edabit = new EdaBit(n);

  // Create local
  EdaBit* edabit0 = new EdaBit(n);
  EdaBit* edabit1 = new EdaBit(n);
  makeLocalEdaBit(edabit0, edabit1, n);

  // Exchange
  EdaBit* tmp_edabit = new EdaBit(n);
  if (server_num == 0) {
    send_EdaBit(serverfd, edabit1, n);
    recv_EdaBit(serverfd, tmp_edabit, n);
    delete edabit1;
    edabit1 = tmp_edabit;
  } else {
    send_EdaBit(serverfd, edabit0, n);
    recv_EdaBit(serverfd, tmp_edabit, n);
    delete edabit0;
    edabit0 = tmp_edabit;
  }
  // Add arithmetic shares
  fmpz_add(edabit->r, edabit0->r, edabit1->r);
  fmpz_mod(edabit->r, edabit->r, Int_Modulus);

  // Add binary shares via circuit
  bool carry = addBinaryShares(serverfd, server_num, n, edabit0->b, edabit1->b, edabit->b, triples);

  // Convert carry to arithmetic [b_n]_p
  fmpz_t bpn;
  fmpz_init(bpn);
  b2a_daBit(serverfd, server_num, dabit, carry, bpn);

  // Subtract out 2^n * [b_n]_p from r
  fmpz_submul_ui(edabit->r, bpn, 1ul << n);
  fmpz_mod(edabit->r, edabit->r, Int_Modulus);

  fmpz_clear(bpn);
  delete edabit0;
  delete edabit1;

  return edabit;
}