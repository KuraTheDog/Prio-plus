#ifndef HE_TRIPLES_H
#define HE_TRIPLES_H

#include <iostream>
#include <random>
#include <vector>

#include "palisade.h"

#include "constants.h"
#include "share.h"

using namespace lbcrypto;

/*
A generator for Arithmetic triples, between two servers.

Uses PALISADE's BFVrns encryption scheme.
The plaintextModulus is the PARSEC Int_Modulus.
  Only supports < 60 bits.
  Due to noise, has only been tested for up to ~48 bits

Can generate up to 8192 triples at a time.
*/
class ArithTripleGenerator {
private:
  const slong plaintextModulus = fmpz_get_si(Int_Modulus);
  const double sigma = 3.2;
  const SecurityLevel securityLevel = HEStd_128_classic;

  CryptoContext<DCRTPoly> cc;

  LPPublicKey<DCRTPoly> pk;
  LPPrivateKey<DCRTPoly> sk;
  LPPublicKey<DCRTPoly> other_pk;

  const int serverfd;

  std::default_random_engine generator;
  std::function<int64_t()> random_int;

  void swapPK();

public:

  /*
   - serverfd: sockfd of the other server
   - server_num: index of the server. Only used for randomness offset
   - random_offset: Have server 1 make this many extra random values, so that the servers have different random values even when starting with the same seed
  */
  ArithTripleGenerator(const int serverfd, const int server_num = 0, const unsigned int random_offset = 8)
  : serverfd(serverfd)
  {
    if (fmpz_cmp_ui(Int_Modulus, 1ULL << 60) > 0) {
      perror("ERROR: PALISADE based triples don't support Int_Modulus >60 bits");
      return;
    } else if (fmpz_cmp_ui(Int_Modulus, 1ULL << 49) > 0) {
      std::cout << "WARNING: PALISADE based triples is not fully tested on Int_Modulus > 48 bits. Running anyways." << std::endl;
    }

    cc = CryptoContextFactory<DCRTPoly>::genCryptoContextBFVrns(
      plaintextModulus, securityLevel, sigma, 0, 1, 0, OPTIMIZED);
    cc->Enable(ENCRYPTION);
    cc->Enable(SHE);
    LPKeyPair<DCRTPoly> keyPair = cc->KeyGen();
    pk = keyPair.publicKey;
    sk = keyPair.secretKey;
    cc->EvalMultKeyGen(sk);

    // TODO: Use FMPZ random instead?
    std::uniform_int_distribution<int64_t> index_dist{0, plaintextModulus - 1};
    random_int = std::bind(index_dist, generator);

    // Randomness offset
    if (server_num == 1) {
      for (unsigned int i = 0; i < random_offset; i++) {
        random_int();
      }
    }

    swapPK();
  }

  std::vector<BeaverTriple*> generateTriples(const size_t n) const;
};

#endif
