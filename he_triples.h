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

  // TODO: Use FMPZ random instead?
  std::default_random_engine generator;
  std::function<int64_t()> random_int;

public:

  /*
   - serverfd: sockfd of the other server
   - server_num: index of the server. Only used for randomness offset
   - random_offset: Have server 1 make this many extra random values, so that the servers have different random values even when starting with the same seed
  */
  ArithTripleGenerator(const int serverfd, const int server_num = 0, const unsigned int random_offset = 8);

  // Make n (up to 8192) arithmetic beaver triples at once
  std::vector<BeaverTriple*> generateTriples(const size_t n) const;
};

#endif
