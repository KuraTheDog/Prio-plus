#ifndef HE_TRIPLES_H
#define HE_TRIPLES_H

#include <iostream>
#include <random>
#include <vector>

#include "palisade.h"

#include "constants.h"
#include "share.h"

using namespace lbcrypto;

// Can make fit at most 8192 entries in a plaintext
#define MAX_HE_BATCH 8192

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

  // Random number generation
  std::default_random_engine generator;
  std::function<int64_t()> random_int;

  // True to work correctly on large batches, and faster
  // False to do debugging
  const bool do_fork;

  // Sends out T array mine, returns other's similar T array
  // T is some serializable object, so public key, cipher text, etc.
  template <class T>
  T* serializedSwap(const size_t num_batches, const T* mine) const;

public:
  ArithTripleGenerator(const int serverfd, const bool do_fork = true);

  // Make n arithmetic beaver triples at once, in batches of 8192
  std::vector<BeaverTriple*> generateTriples(const size_t n) const;
};

#endif
