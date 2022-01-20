#include "he_triples.h"

#include <sys/wait.h>

#include <string>

// Serializing
#include "ciphertext-ser.h"
#include "pubkeylp-ser.h"
#include "scheme/bfvrns/bfvrns-ser.h"
#include "utils/serial.h"

#include "net_share.h"
#include "share.h"

// Sends a serializable T mine, receives sent T other.
template <class T>
T* ArithTripleGenerator::serializedSwap(const size_t num_batches, const T* mine) const {
  pid_t pid = 0;
  int status = 0;
  if (do_fork) pid = fork();
  if (pid == 0) {
    for (unsigned int i = 0; i < num_batches; i++) {
      std::string s;
      std::stringstream ss;
      Serial::Serialize(mine[i], ss, SerType::BINARY);
      s = ss.str();
      send_string(serverfd, s);
    }

    if (do_fork) exit(EXIT_SUCCESS);
  }

  T* other = new T[num_batches];
  for (unsigned int i = 0; i < num_batches; i++) {
    std::string s2;
    std::stringstream ss2;
    recv_string(serverfd, s2);
    ss2 << s2;
    ss2.flush();
    Serial::Deserialize(other[i], ss2, SerType::BINARY);
  }

  if (do_fork) waitpid(pid, &status, 0);
  return other;
}

ArithTripleGenerator::ArithTripleGenerator(const int serverfd,
                                           const bool do_fork)
: serverfd(serverfd)
, do_fork(do_fork)
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

  // Seed and set up random ints
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int64_t> index_dist{0, plaintextModulus - 1};
  random_int = std::bind(index_dist, generator);

  // Swap public keys
  LPPublicKey<DCRTPoly> inp[1] = {pk};
  LPPublicKey<DCRTPoly>* tmp = serializedSwap(1, inp);
  other_pk = tmp[0];
  delete[] tmp;
}

std::vector<BeaverTriple*> ArithTripleGenerator::generateTriples(const size_t n) const {
  std::vector<BeaverTriple*> res;
  if (n == 0)
    return res;
  res.reserve(n);

  // ceil(n / MAX_HE_BATCH)
  const size_t num_batches = (n - 1) / MAX_HE_BATCH + 1;
  const size_t last_batch = (n - 1) % MAX_HE_BATCH + 1;

  std::vector<int64_t>* a = new std::vector<int64_t>[num_batches];
  std::vector<int64_t>* b = new std::vector<int64_t>[num_batches];
  std::vector<int64_t>* d = new std::vector<int64_t>[num_batches];

  Ciphertext<DCRTPoly>* ct_a = new Ciphertext<DCRTPoly>[num_batches];
  Ciphertext<DCRTPoly>* ct_d = new Ciphertext<DCRTPoly>[num_batches];
  Plaintext* plain_b = new Plaintext[num_batches];

  for (unsigned int i = 0; i < num_batches; i++) {
    const size_t N = (i == num_batches - 1 ? last_batch : MAX_HE_BATCH);
    // local random values
    a[i].reserve(N);
    b[i].reserve(N);
    d[i].reserve(N);
    for (unsigned int j = 0; j < N; j++) {
      a[i].push_back(random_int());
      b[i].push_back(random_int());
      d[i].push_back(random_int());
    }

    Plaintext plain_a = cc->MakePackedPlaintext(a[i]);
    plain_b[i] = cc->MakePackedPlaintext(b[i]);
    Plaintext plain_d = cc->MakePackedPlaintext(d[i]);
    ct_a[i] = cc->Encrypt(pk, plain_a);        // Enc(a)
    ct_d[i] = cc->Encrypt(other_pk, plain_d);  // Enc'(d)
  }

  // Swap Enc(a), get Enc'(a')
  Ciphertext<DCRTPoly>* ct_a2 = serializedSwap(num_batches, ct_a);

  Ciphertext<DCRTPoly>* ct_e2 = new Ciphertext<DCRTPoly>[num_batches];
  for (unsigned int i = 0; i < num_batches; i++) {
    // E' = b Enc'(a') - Enc'(d)
    auto ct_a2_b = cc->EvalMult(ct_a2[i], plain_b[i]);
    ct_e2[i] = cc->EvalSub(ct_a2_b, ct_d[i]);
  }

  // Swap E, get E = b' Enc(a) - Enc(d')
  Ciphertext<DCRTPoly>* ct_e = serializedSwap(num_batches, ct_e2);

  for (unsigned int i = 0; i < num_batches; i++) {
    const size_t N = (i == num_batches - 1 ? last_batch : MAX_HE_BATCH);
    // e = Dec(E) = b' a - d', so d' + e = a b'
    Plaintext plain_e;
    cc->Decrypt(sk, ct_e[i], &plain_e);
    // Truncate e to (relevant) first N values;
    plain_e->SetLength(N);
    auto e = plain_e->GetPackedValue();

    // Build result
    for (unsigned int j = 0; j < N; j++) {
      BeaverTriple* triple = new BeaverTriple();
      fmpz_set_si(triple->A, a[i][j]);
      fmpz_set_si(triple->B, b[i][j]);
      // c = ab + d + e
      fmpz_set_si(triple->C, a[i][j]);
      fmpz_mul_si(triple->C, triple->C, b[i][j]);
      fmpz_mod(triple->C, triple->C, Int_Modulus);
      fmpz_add_si(triple->C, triple->C, d[i][j]);
      fmpz_add_si(triple->C, triple->C, e[j]);
      fmpz_mod(triple->C, triple->C, Int_Modulus);

      res.push_back(triple);
    }
  }

  delete[] a;
  delete[] b;
  delete[] d;
  delete[] plain_b;
  delete[] ct_a;
  delete[] ct_d;
  delete[] ct_a2;
  delete[] ct_e2;
  delete[] ct_e;

  return res;
}
