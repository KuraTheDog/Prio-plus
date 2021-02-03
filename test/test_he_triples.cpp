#include <iostream>
#include <random>
#include <string>

#include "palisade.h"

#include "utils_test_connect.h"
#include "../share.h"
#include "../net_share.h"
#include "../he_triples.h"
// #include "../constants.h"

using namespace lbcrypto;

// Basic run: local, one context/keys. 
void runLocal() {

  int N = 20;

  int64_t plaintextModulus = 65537;  // q = 1 mod 2^14 = m
  double sigma = 3.2;
  SecurityLevel securityLevel = HEStd_128_classic;
  uint32_t depth = 1;

  CryptoContext<DCRTPoly> cryptoContext =
    CryptoContextFactory<DCRTPoly>::genCryptoContextBFVrns(
      plaintextModulus, securityLevel, sigma, 0, depth, 0, OPTIMIZED);

  cryptoContext->Enable(ENCRYPTION);
  cryptoContext->Enable(SHE);
  // cryptoContext->Enable(LEVELEDSHE);

  std::default_random_engine generator;
  // Pick a random coefficient index from 0 to N-1
  std::uniform_int_distribution<int64_t> index_dist{0, plaintextModulus};
  std::function<int64_t()> rint = std::bind(index_dist, generator);

  // Initialize Public Key Containers
  LPKeyPair<DCRTPoly> keyPair;

  // Generate a public/private key pair
  keyPair = cryptoContext->KeyGen();

  // Generate the relinearization key
  cryptoContext->EvalMultKeyGen(keyPair.secretKey);

  std::vector<int64_t> a;
  std::vector<int64_t> b;
  std::vector<int64_t> d;
  for (int i = 0; i < N; i++) {
    a.push_back(rint()); b.push_back(rint()); d.push_back(rint());
  }
  Plaintext plain_a = cryptoContext->MakePackedPlaintext(a);
  Plaintext plain_b = cryptoContext->MakePackedPlaintext(b);
  Plaintext plain_d = cryptoContext->MakePackedPlaintext(d);
  Ciphertext<DCRTPoly> ct_a = cryptoContext->Encrypt(keyPair.publicKey, plain_a);
  Ciphertext<DCRTPoly> ct_d = cryptoContext->Encrypt(keyPair.publicKey, plain_d);

  std::vector<int64_t> a2;
  std::vector<int64_t> b2;
  std::vector<int64_t> d2;
  for (int i = 0; i < N; i++) {
    a2.push_back(rint()); b2.push_back(rint()); d2.push_back(rint());
  }
  Plaintext plain_a2 = cryptoContext->MakePackedPlaintext(a2);
  Plaintext plain_b2 = cryptoContext->MakePackedPlaintext(b2);
  Plaintext plain_d2 = cryptoContext->MakePackedPlaintext(d2);
  auto ct_a2 = cryptoContext->Encrypt(keyPair.publicKey, plain_a2);
  auto ct_d2 = cryptoContext->Encrypt(keyPair.publicKey, plain_d2);

  // Swap a -> 2, a2 -> 1

  // Server 1 calc
  auto ct_a2_b = cryptoContext->EvalMult(ct_a2, plain_b);
  auto ct_e2 = cryptoContext->EvalAdd(ct_a2_b, ct_d);

  // Server 2 calc
  auto ct_a_b2 = cryptoContext->EvalMult(ct_a, plain_b2);
  auto ct_e = cryptoContext->EvalAdd(ct_a_b2, ct_d2);  

  // swap e2 -> 2, e -> 1

  // Decrypt the result of additions
  Plaintext plain_e;
  cryptoContext->Decrypt(keyPair.secretKey, ct_e, &plain_e);
  plain_e->SetLength(N);
  auto e = plain_e->GetPackedValue();

  Plaintext plain_e2;
  cryptoContext->Decrypt(keyPair.secretKey, ct_e2, &plain_e2);
  plain_e2->SetLength(N);
  auto e2 = plain_e2->GetPackedValue();

  std::cout << "a : " << a << std::endl;
  std::cout << "a2: " << a2 << std::endl;
  std::cout << "b : " << b << std::endl;
  std::cout << "b2: " << b2 << std::endl;

  std::cout << std::endl;

  std::cout << "d : " << d << std::endl;
  std::cout << "d2: " << d2 << std::endl;

  std::cout << std::endl;  

  std::cout << "e : " << e << std::endl;
  std::cout << "e2: " << e2 << std::endl;

  for (int i = 0; i < N; i++) {
    std::cout << i << ": ";
    int64_t c = (a[i] * b[i] - d[i] + e[i]) % plaintextModulus;
    int64_t c2 = (a2[i] * b2[i] - d2[i] + e2[i]) % plaintextModulus;
    std::cout << "(" << a[i] << " + " << a2[i] << ")";
    std::cout << " * (" << b[i] << " + " << b2[i] << ")";
    std::cout << " = " << ((a[i] + a2[i]) * (b[i] + b2[i])) % plaintextModulus;
    std::cout << ", (" << c << " + " << c2 << ")";
    std::cout << " = " << (c + c2) % plaintextModulus << std::endl;
  }
}


void runServerTest(const int server_num, const int serverfd, const size_t N) {
  auto start = emp::clock_start();

  // int64_t modulus = 65537;
  int64_t modulus = 0x8008001;
  // int64_t modulus = 0x800008001L;
  // int64_t modulus = 0x80000000080001L;

  ArithTripleGenerator gen(serverfd, server_num, modulus);
  std::cout << server_num << " initialized in: " << (((float) emp::time_from(start)) / CLOCKS_PER_SEC) << "s" << std::endl; start = emp::clock_start();
  gen.swapPK();
  std::cout << server_num << " PK Swapped in: " << (((float) emp::time_from(start)) / CLOCKS_PER_SEC) << "s" << std::endl; start = emp::clock_start();

  std::vector<BeaverTriple*> triples = gen.generateTriples(N);

  std::cout << server_num << " made # triples: " << triples.size() << std::endl;

  std::cout << server_num << " generated in: " << (((float) emp::time_from(start)) / CLOCKS_PER_SEC) << "s" << std::endl;

  std::cout << server_num << " triple 0: ";
  fmpz_print(triples[0]->A); std::cout << ", ";
  fmpz_print(triples[0]->B); std::cout << ", ";
  fmpz_print(triples[0]->C); std::cout << std::endl;

  std::cout << server_num << " triple " << N-1 << ": ";
  fmpz_print(triples[N-1]->A); std::cout << ", ";
  fmpz_print(triples[N-1]->B); std::cout << ", ";
  fmpz_print(triples[N-1]->C); std::cout << std::endl;

  fmpz_t mod;
  fmpz_init_set_si(mod, modulus);
  std::cout << "modulus: ";
  fmpz_print(mod);
  std::cout << std::endl;

  int idx = 0;
  std::cout << "Validating triple: " << idx << std::endl;
  if (server_num == 0) {
    BeaverTriple* triple = triples[idx];
    BeaverTriple* other_triple = new BeaverTriple();
    recv_BeaverTriple(serverfd, other_triple);

    fmpz_t tmp; fmpz_init(tmp);
    fmpz_t tmp2; fmpz_init(tmp2);
    fmpz_add(tmp, triple->A, other_triple->A);
    fmpz_mod(tmp, tmp, mod);
    fmpz_add(tmp2, triple->B, other_triple->B);
    fmpz_mod(tmp2, tmp2, mod);
    fmpz_mul(tmp, tmp, tmp2);
    fmpz_mod(tmp, tmp, mod);
    std::cout << "actual product: ("; fmpz_print(triple->A);
    std::cout << " + "; fmpz_print(other_triple->A);
    std::cout << ") * ("; fmpz_print(triple->B);
    std::cout << " + "; fmpz_print(other_triple->B);
    std::cout << ") = "; fmpz_print(tmp);
    std::cout << std::endl;
    fmpz_add(tmp, triple->C, other_triple->C);
    fmpz_mod(tmp, tmp, mod);
    fmpz_print(triple->C); std::cout << " + "; fmpz_print(other_triple->C);
    std::cout << " = "; fmpz_print(tmp); std::cout << std::endl;
  } else {
    send_BeaverTriple(serverfd, triples[idx]);
    // send_BooleanBeaverTriple(newsockfd, triples[N-1]);
  }
}

void serverTest(const size_t N) {
  int sockfd = init_receiver();

  if (fork() == 0) {
    int cli_sockfd = init_sender();
    runServerTest(0, cli_sockfd, N);
    close(cli_sockfd);
  } else {
    int newsockfd = accept_receiver(sockfd);
    runServerTest(1, newsockfd, N);
    close(newsockfd);
  }
  close(sockfd);
}


int main(int argc, char** argv){
  // runLocal();
  serverTest(10);

  return 0;
}