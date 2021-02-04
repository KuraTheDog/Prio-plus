#include "he_triples.h"

#include <string>

// Serializing
#include "ciphertext-ser.h"
#include "pubkeylp-ser.h"
#include "scheme/bfvrns/bfvrns-ser.h"
#include "utils/serialize-binary.h"

#include "net_share.h"
#include "share.h"

// Sends a serializable T mine, recieves sent T other.
template <class T>
T serializedSwap(const int serverfd, const T mine) {
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

void ArithTripleGenerator::swapPK() {
  other_pk = serializedSwap(serverfd, pk);
}

std::vector<BeaverTriple*> ArithTripleGenerator::generateTriples(const size_t n) const {
  if (n > 8192) {
    std::cout << "Can only make 8192 triples at a time" << std::endl;
  }
  const size_t N = (n > 8192 ? 8192 : n);

  // Generate random local values
  std::vector<int64_t> a, b, d;
  for (int i = 0; i < N; i++) {
    a.push_back(random_int());
    b.push_back(random_int());
    d.push_back(random_int());
  }

  Plaintext plain_a = cc->MakePackedPlaintext(a);
  Plaintext plain_b = cc->MakePackedPlaintext(b);
  Plaintext plain_d = cc->MakePackedPlaintext(d);
  auto ct_a = cc->Encrypt(pk, plain_a);          // Enc(a)
  auto ct_d = cc->Encrypt(other_pk, plain_d);    // Enc'(d)

  // Swap a: get Enc'(a')
  Ciphertext<DCRTPoly> ct_a2 = serializedSwap(serverfd, ct_a);

  // E' = b Enc'(a') - Enc'(d)
  auto ct_a2_b = cc->EvalMult(ct_a2, plain_b);
  auto ct_e2 = cc->EvalSub(ct_a2_b, ct_d);

  // Swap E, get E = b' Enc(a) - Enc(d')
  Ciphertext<DCRTPoly> ct_e = serializedSwap(serverfd, ct_e2);

  // e = Dec(E) = b' a - d', so d' + e = a b'
  Plaintext plain_e;
  cc->Decrypt(sk, ct_e, &plain_e);
  // Truncate e to (relevant) first n values
  plain_e->SetLength(N);
  auto e = plain_e->GetPackedValue();

  std::vector<BeaverTriple*> res;
  for (unsigned int i = 0; i < N; i++) {
    BeaverTriple* triple = new BeaverTriple();
    fmpz_set_si(triple->A, a[i]);
    fmpz_set_si(triple->B, b[i]);
    // c = ab + d + e
    fmpz_set_si(triple->C, a[i]);
    fmpz_mul_si(triple->C, triple->C, b[i]);
    fmpz_add_si(triple->C, triple->C, d[i]);
    fmpz_add_si(triple->C, triple->C, e[i]);
    fmpz_mod(triple->C, triple->C, Int_Modulus);

    res.push_back(triple);
  }

  return res;
}
