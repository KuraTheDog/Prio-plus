#include <iostream>
#include <random>

#include "palisade.h"

#include "../he_triples.h"

#include "utils_test_connect.h"
#include "../constants.h"
#include "../net_share.h"
#include "../share.h"
#include "../utils.h"

using namespace lbcrypto;

// Basic run: local, one context/keys
void runLocal(size_t N) {

  // Does account for max batch size
  if (N > 8192) N = 8192;

  // int64_t plaintextModulus = 65537;
  // int64_t plaintextModulus = 0x8008001;
  int64_t plaintextModulus = 0x800008001;
  double sigma = 3.2;
  SecurityLevel securityLevel = HEStd_128_classic;
  uint32_t depth = 1;

  // Must be between 30 and 60.
  // must be greater than Plaintext modulus
  // For 0x800008001:
  // 35 is too small. 36 doesn't work
  unsigned int dcrtBits = 60;

  CryptoContext<DCRTPoly> cryptoContext =
    CryptoContextFactory<DCRTPoly>::genCryptoContextBFVrns(
      plaintextModulus,
      securityLevel,
      sigma,
      0,
      depth,
      0,
      OPTIMIZED,
      2,
      0,
      dcrtBits
  );

  cryptoContext->Enable(ENCRYPTION);
  cryptoContext->Enable(SHE);
  // cryptoContext->Enable(LEVELEDSHE);

  int64_t maxrand = 1L << 50;
  maxrand = (plaintextModulus < maxrand ? plaintextModulus : maxrand) - 1;

  std::cout << "plaintextModulus: " << plaintextModulus << std::endl;

  std::default_random_engine generator;
  // Pick a random coefficient index from 0 to N-1
  std::uniform_int_distribution<int64_t> index_dist{0, maxrand};
  std::function<int64_t()> rint = std::bind(index_dist, generator);

  // Initialize Public Key Containers
  LPKeyPair<DCRTPoly> keyPair;

  // Generate a public/private key pair
  keyPair = cryptoContext->KeyGen();

  // Generate the relinearization key
  cryptoContext->EvalMultKeyGen(keyPair.secretKey);
  // cryptoContext->EvalSumKeyGen(keyPair.secretKey);

  std::vector<int64_t> a;
  std::vector<int64_t> b;
  std::vector<int64_t> d;
  for (unsigned int i = 0; i < N; i++) {
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
  for (unsigned int i = 0; i < N; i++) {
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
  auto ct_e2 = cryptoContext->EvalSub(ct_a2_b, ct_d);

  // Server 2 calc
  auto ct_a_b2 = cryptoContext->EvalMult(ct_a, plain_b2);
  auto ct_e = cryptoContext->EvalSub(ct_a_b2, ct_d2);

  Plaintext plain_tmp;
  cryptoContext->Decrypt(keyPair.secretKey, ct_a2_b, &plain_tmp);
  plain_tmp->SetLength(N); auto tmp_val = plain_tmp->GetPackedValue();
  // std::cout << "dec(a2_b): " << tmp_val << std::endl;
  cryptoContext->Decrypt(keyPair.secretKey, ct_e2, &plain_tmp);
  plain_tmp->SetLength(N); tmp_val = plain_tmp->GetPackedValue();
  // std::cout << "dec(e2):   " << tmp_val << std::endl;

  cryptoContext->Decrypt(keyPair.secretKey, ct_a_b2, &plain_tmp);
  plain_tmp->SetLength(N); tmp_val = plain_tmp->GetPackedValue();
  // std::cout << "dec(a_b2): " << tmp_val << std::endl;
  cryptoContext->Decrypt(keyPair.secretKey, ct_e, &plain_tmp);
  plain_tmp->SetLength(N); tmp_val = plain_tmp->GetPackedValue();
  // std::cout << "dec(e):    " << tmp_val << std::endl;

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

  // std::cout << "a : " << a << std::endl;
  // std::cout << "a2: " << a2 << std::endl;
  // std::cout << "b : " << b << std::endl;
  // std::cout << "b2: " << b2 << std::endl;

  // std::cout << std::endl;

  // std::cout << "d : " << d << std::endl;
  // std::cout << "d2: " << d2 << std::endl;

  // std::cout << std::endl;

  // std::cout << "e : " << e << std::endl;
  // std::cout << "e2: " << e2 << std::endl;

  fmpz_t c_val; fmpz_init(c_val);
  fmpz_t c2_val; fmpz_init(c2_val);
  fmpz_t fmod; fmpz_init_set_si(fmod, plaintextModulus);
  fmpz_t tmp; fmpz_init(tmp);
  fmpz_t tmp2; fmpz_init(tmp2);

  for (unsigned int i = 0; i < N; i++) {
    // std::cout << i << ": ";
    fmpz_set_si(c_val, a[i]); fmpz_mul_si(c_val, c_val, b[i]);
    fmpz_add_si(c_val, c_val, d[i]); fmpz_add_si(c_val, c_val, e[i]);
    fmpz_mod(c_val, c_val, fmod);

    fmpz_set_si(c2_val, a2[i]); fmpz_mul_si(c2_val, c2_val, b2[i]);
    fmpz_add_si(c2_val, c2_val, d2[i]); fmpz_add_si(c2_val, c2_val, e2[i]);
    fmpz_mod(c2_val, c2_val, fmod);

    fmpz_set_si(tmp, a[i]); fmpz_add_si(tmp, tmp, a2[i]);
    fmpz_set_si(tmp2, b[i]); fmpz_add_si(tmp2, tmp2, b2[i]);
    fmpz_mul(tmp, tmp, tmp2); fmpz_mod(tmp, tmp, fmod);
    fmpz_add(tmp2, c_val, c2_val); fmpz_mod(tmp2, tmp2, fmod);
    bool match = fmpz_equal(tmp, tmp2);
    if (i == 0 and !match) {
      std::cout << i << " does not match" << std::endl;
      std::cout << "[a] * [b] = (" << a[i] << " + " << a2[i] << ")";
      std::cout << " * (" << b[i] << " + " << b2[i] << ")";
      std::cout << "\n = "; fmpz_print(tmp);
      std::cout << "\n  [c] = ("; fmpz_print(c_val);
      std::cout << " + "; fmpz_print(c_val);
      std::cout << ")\n = "; fmpz_print(tmp2); std::cout << std::endl;
    }
  }
}

void runServerTest(const int server_num, const int serverfd, const size_t N) {
  auto start = clock_start();

  // int64_t modulus = 0x8008001;           // 26 bit
  // int64_t modulus = 0x800008001L;           // 36 bit, works up to 1<<32
  // int64_t modulus = 0x80000000080001LL;  // 55 bit, doesn't run
  // int64_t modulus = 0x8000002080001LL;
  // int64_t modulus = 0x800006880001LL;
  // std::cout << "modulus:    " << modulus << std::endl;

  ArithTripleGenerator gen(serverfd, server_num);
  std::cout << server_num << " initialized gen in: " << sec_from(start) << "s" << std::endl; start = clock_start();

  std::vector<BeaverTriple*> triples = gen.generateTriples(N);

  std::cout << server_num << " generated " << triples.size() << " triples in: ";
  std::cout << sec_from(start) << "s" << std::endl; start = clock_start();

  std::cout << server_num << "'s triple 0: ";
  fmpz_print(triples[0]->A); std::cout << ", ";
  fmpz_print(triples[0]->B); std::cout << ", ";
  fmpz_print(triples[0]->C); std::cout << std::endl;

  if (N > 1) {
    std::cout << server_num << "'s triple " << N-1 << ": ";
    fmpz_print(triples[N-1]->A); std::cout << ", ";
    fmpz_print(triples[N-1]->B); std::cout << ", ";
    fmpz_print(triples[N-1]->C); std::cout << std::endl;
  }

  for (unsigned int idx = 0; idx < N; idx++) {
    // std::cout << "Validating triple: " << idx << std::endl;
    if (server_num == 0) {
      BeaverTriple* triple = triples[idx];
      BeaverTriple* other_triple = new BeaverTriple();
      recv_BeaverTriple(serverfd, other_triple);

      fmpz_t tmp; fmpz_init(tmp);
      fmpz_t tmp2; fmpz_init(tmp2);
      fmpz_add(tmp, triple->A, other_triple->A);
      fmpz_mod(tmp, tmp, Int_Modulus);
      fmpz_add(tmp2, triple->B, other_triple->B);
      fmpz_mod(tmp2, tmp2, Int_Modulus);
      fmpz_mul(tmp, tmp, tmp2);
      fmpz_mod(tmp, tmp, Int_Modulus);
      fmpz_add(tmp2, triple->C, other_triple->C);
      fmpz_mod(tmp2, tmp2, Int_Modulus);
      bool valid = fmpz_equal(tmp, tmp2);
      send_bool(serverfd, valid);
      if (not valid or idx == N-1) {
        if (not valid) {
          std::cout << "############## Invalid " << idx << " ##############" << std::endl;
        } else {
          std::cout << "All valid. Math of N-1: " << std::endl;
        }
        std::cout << "actual product: ("; fmpz_print(triple->A);
        std::cout << " + "; fmpz_print(other_triple->A);
        std::cout << ") * ("; fmpz_print(triple->B);
        std::cout << " + "; fmpz_print(other_triple->B);
        std::cout << ") = \n"; fmpz_print(tmp);
        std::cout << std::endl;
        fmpz_print(triple->C); std::cout << " + "; fmpz_print(other_triple->C);
        std::cout << " = \n"; fmpz_print(tmp2); std::cout << std::endl;
      }
      delete other_triple;
      delete triples[idx];
      if (not valid) break;
    } else {
      send_BeaverTriple(serverfd, triples[idx]);
      bool valid;
      recv_bool(serverfd, valid);
      delete triples[idx];
      if (not valid) break;
    }
  }

  std::cout << server_num << " ran all validation in: " << sec_from(start) << "s" << std::endl;
}

void serverTest(const size_t N) {

  std::thread t0([&]() {
    int cli_sockfd = init_sender();
    runServerTest(0, cli_sockfd, N);
    close(cli_sockfd);
  });
  std::thread t1([&]() {
    int sockfd = init_receiver();
    int newsockfd = accept_receiver(sockfd);
    runServerTest(1, newsockfd, N);
    close(newsockfd);
    close(sockfd);
  });
  t0.join();
  t1.join();
}


int main(int argc, char** argv){
  init_constants();

  int server_num = -1;
  if(argc >= 2){
    server_num = atoi(argv[1]);
  }

  const size_t N = 20000;

  // runLocal(N);

  if (server_num == -1) {
    serverTest(N);
  } else if (server_num == 0) {
    int sockfd = init_receiver();
    int newsockfd = accept_receiver(sockfd);
    runServerTest(0, newsockfd, N);
    close(newsockfd);
    close(sockfd);
  } else if (server_num == 1) {
    int cli_sockfd = init_sender();
    runServerTest(1, cli_sockfd, N);
    close(cli_sockfd);
  }


  clear_constants();
  return 0;
}
