#include "he_triples.h"

#include <string>

// Serializing
#include "ciphertext-ser.h"
#include "pubkeylp-ser.h"
#include "scheme/bfvrns/bfvrns-ser.h"
#include "utils/serialize-binary.h"

#include "net_share.h"
#include "share.h"

/*
template class<T>;
T swap(const T mine) {
  std::string s;
  std::stringstream ss;
  Serial::Serialize(mine, ss, SerType::BINARY);
  s = ss.str();
  send_string(serverfd, s);

  T other;
  std::string s2;
  std::stringstream ss2;
  recv_string(serverfd, s2);
  ss2 << s2;
  ss2.flush();
  Serial::Deserialize(other, ss2, SerType::BINARY);

  return other;
}
*/

void ArithTripleGenerator::swapPK() {
  std::string s;
  std::stringstream ss;
  Serial::Serialize(pk, ss, SerType::BINARY);
  s = ss.str();
  send_string(serverfd, s);

  std::string s2;
  std::stringstream ss2;
  recv_string(serverfd, s2);
  ss2 << s2;
  ss2.flush();
  Serial::Deserialize(other_pk, ss2, SerType::BINARY);
}

Ciphertext<DCRTPoly> ArithTripleGenerator::swapCT(const Ciphertext<DCRTPoly> ct) {
  std::string s;
  std::stringstream ss;
  Serial::Serialize(ct, ss, SerType::BINARY);
  s = ss.str();
  int n = send_string(serverfd, s);
  if (debug) {
    std::cout << "sent ct in bits: " << n <<  std::endl;
    debug = false;
  }

  Ciphertext<DCRTPoly> ct2;
  std::string s2;
  std::stringstream ss2;
  recv_string(serverfd, s2);
  ss2 << s2;
  ss2.flush();
  Serial::Deserialize(ct2, ss2, SerType::BINARY);
  return ct2;
}

std::vector<BeaverTriple*> ArithTripleGenerator::generateTriples(const size_t n) {
  // Generate random local values
  std::vector<int64_t> a, b, d;
  for (int i = 0; i < n; i++) {
    a.push_back(random_int());
    b.push_back(random_int());
    d.push_back(random_int());
  }

  Plaintext plain_a = cc->MakePackedPlaintext(a);
  Plaintext plain_b = cc->MakePackedPlaintext(b);
  Plaintext plain_d = cc->MakePackedPlaintext(d);
  auto ct_a = cc->Encrypt(pk, plain_a);  // Enc(a)
  auto ct_d = cc->Encrypt(other_pk, plain_d);    // Enc'(d)

  // Swap a: get Enc'(a')
  Ciphertext<DCRTPoly> ct_a2 = swapCT(ct_a);

  // E' = b Enc'(a') - Enc'(d)
  auto ct_a2_b = cc->EvalMult(ct_a2, plain_b);
  auto ct_e2 = cc->EvalSub(ct_a2_b, ct_d);

  // Swap E, get E = b' Enc(a) - Enc(d')
  Ciphertext<DCRTPoly> ct_e = swapCT(ct_e2);

  // e = Dec(E) = b' a - d', so d' + e = a b'
  Plaintext plain_e;
  cc->Decrypt(sk, ct_e, &plain_e);
  // Truncate e to (relevant) first n values
  plain_e->SetLength(n);
  auto e = plain_e->GetPackedValue();

  std::vector<BeaverTriple*> res;
  for (int i = 0; i < n; i++) {
    BeaverTriple* triple = new BeaverTriple();
    fmpz_set_si(triple->A, a[i]);
    fmpz_set_si(triple->B, b[i]);
    // c = ab + d + e
    fmpz_set_si(triple->C, (a[i] * b[i] + d[i] + e[i]) % plaintextModulus );

    res.push_back(triple);
  }

  return res;
}
