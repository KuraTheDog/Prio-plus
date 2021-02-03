#ifndef HE_TRIPLES_H
#define HE_TRIPLES_H
#include <iostream>
#include <random>
#include <vector>

#include "palisade.h"

#include "share.h"

using namespace lbcrypto;

class ArithTripleGenerator {
  // const int64_t plaintextModulus = 34359771137;  // TODO: Int_Modulus.
  const int64_t plaintextModulus;  // TODO: Int_Modulus.
  const double sigma = 3.2;
  const SecurityLevel securityLevel = HEStd_128_classic;

  CryptoContext<DCRTPoly> cc;

  LPPublicKey<DCRTPoly> pk;
  LPPrivateKey<DCRTPoly> sk;
  LPPublicKey<DCRTPoly> other_pk;

  const int serverfd;

  std::default_random_engine generator;
  std::function<int64_t()> random_int;

  Ciphertext<DCRTPoly> swapCT(const Ciphertext<DCRTPoly> ct);

  bool debug = true;

public:

  ArithTripleGenerator(const int serverfd, const int server_num = 0, const int64_t mod = 65537) 
  : plaintextModulus(mod)
  , serverfd(serverfd) 
  {
    cc = CryptoContextFactory<DCRTPoly>::genCryptoContextBFVrns(
      plaintextModulus, securityLevel, sigma, 0, 1, 0, OPTIMIZED);
    cc->Enable(ENCRYPTION);
    cc->Enable(SHE);
    LPKeyPair<DCRTPoly> keyPair = cc->KeyGen();
    pk = keyPair.publicKey;
    sk = keyPair.secretKey;
    cc->EvalMultKeyGen(sk);

    std::uniform_int_distribution<int64_t> index_dist{0, plaintextModulus};
    random_int = std::bind(index_dist, generator);

    // Randomness offset
    if (server_num == 1) {
      random_int();
    }
  }

  void swapPK();

  std::vector<BeaverTriple*> generateTriples(const size_t n);

  int64_t getMod() const {
    return plaintextModulus;
  }
};

#endif
