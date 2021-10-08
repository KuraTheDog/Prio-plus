#include "correlated.h"

#include <sys/wait.h>

#include <iostream>

#include "constants.h"
#include "net_share.h"
#include "ot.h"
#include "utils.h"

void CorrelatedStore::addBoolTriples(const size_t n) {
  auto start = clock_start();
  const size_t num_to_make = (n > batch_size ? n : batch_size);
  std::cout << "adding booltriples: " << num_to_make << std::endl;
  std::queue<BooleanBeaverTriple*> new_triples = gen_boolean_beaver_triples(server_num, num_to_make, ot0, ot1);
  for (unsigned int i = 0; i < num_to_make; i++) {
    btriple_store.push(new_triples.front());
    new_triples.pop();
  }
  std::cout << "addBoolTriples timing : " << sec_from(start) << std::endl;
}

void CorrelatedStore::addDaBits(const size_t n) {
  auto start = clock_start();
  // const size_t num_to_make = (n > batch_size ? n : batch_size);
  const size_t num_to_make = n;  // Currently to make "end to end" easier to benchmark
  std::cout << "adding dabits: " << num_to_make << std::endl;
  if (!lazy) {
    DaBit** dabit = generateDaBit(num_to_make);
    for (unsigned int i = 0; i < num_to_make; i++)
      dabit_store.push(dabit[i]);
    delete[] dabit;
  } else {  // Lazy generation: make local and send over
    DaBit** const dabit = new DaBit*[num_to_make];
    for (unsigned int i = 0; i < num_to_make; i++)
      dabit[i] = new DaBit;
    if (server_num == 0) {
      DaBit** const other_dabit = new DaBit*[num_to_make];
      for (unsigned int i = 0; i < num_to_make; i++) {
        other_dabit[i] = new DaBit;
        makeLocalDaBit(dabit[i], other_dabit[i]);
      }
      send_DaBit_batch(serverfd, other_dabit, num_to_make);
      for (unsigned int i = 0; i < num_to_make; i++) {
        dabit_store.push(dabit[i]);
        delete other_dabit[i];
      }
      delete[] other_dabit;
    } else {
      recv_DaBit_batch(serverfd, dabit, num_to_make);
      for (unsigned int i = 0; i < num_to_make; i++)
        dabit_store.push(dabit[i]);
    }
    delete[] dabit;
  }
  std::cout << "addDaBits timing : " << sec_from(start) << std::endl;
}

void CorrelatedStore::checkBoolTriples(const size_t n) {
  if (btriple_store.size() < n) addBoolTriples(n - btriple_store.size());
}

void CorrelatedStore::checkDaBits(const size_t n) {
  if (dabit_store.size() < n) addDaBits(n - dabit_store.size());
}

BooleanBeaverTriple* CorrelatedStore::getBoolTriple() {
  checkBoolTriples(1);
  BooleanBeaverTriple* ans = btriple_store.front();
  btriple_store.pop();
  return ans;
}

DaBit* CorrelatedStore::getDaBit() {
  checkDaBits(1);
  DaBit* ans = dabit_store.front();
  dabit_store.pop();
  return ans;
}

void CorrelatedStore::printSizes() {
  std::cout << "Current store sizes:" << std::endl;
  std::cout << " Dabits: " << dabit_store.size() << std::endl;
  // std::cout << " Bool  Triples: " << btriple_store.size() << std::endl;
}

void CorrelatedStore::maybeUpdate() {
  auto start = clock_start();

  // Make top level if stores not enough
  const bool make_da = dabit_store.size() < (batch_size / 2);
  // Determine how much of each to make
  const size_t da_target = 2 * make_da;
  const size_t btrip_target = 0;  // NOTE: Currently disabled

  if (btriple_store.size() < btrip_target * batch_size)
      addBoolTriples(btrip_target * batch_size);

  if (dabit_store.size() < da_target * batch_size)
    addDaBits(da_target * batch_size);

  printSizes();

  std::cout << "precompute timing : " << sec_from(start) << std::endl;
}

CorrelatedStore::~CorrelatedStore() {
  while (!dabit_store.empty()) {
    DaBit* bit = dabit_store.front();
    dabit_store.pop();
    delete bit;
  }
  while (!btriple_store.empty()) {
    BooleanBeaverTriple* triple = btriple_store.front();
    btriple_store.pop();
    delete triple;
  }
}

bool* CorrelatedStore::multiplyBoolShares(const size_t N,
                                          const bool* const x,
                                          const bool* const y) {
  bool* z = new bool[N];

  bool* d_this = new bool[N];
  bool* e_this = new bool[N];

  checkBoolTriples(N);
  for (unsigned int i = 0; i < N; i++) {
    BooleanBeaverTriple* triple = getBoolTriple();
    d_this[i] = x[i] ^ triple->a;
    e_this[i] = y[i] ^ triple->b;
    z[i] = triple->c;
    delete triple;
  }
  pid_t pid = 0;
  int status = 0;
  if (do_fork) pid = fork();
  if (pid == 0) {
    send_bool_batch(serverfd, d_this, N);
    send_bool_batch(serverfd, e_this, N);

    if (do_fork) exit(EXIT_SUCCESS);
  }

  bool* d_other = new bool[N];
  bool* e_other = new bool[N];
  recv_bool_batch(serverfd, d_other, N);
  recv_bool_batch(serverfd, e_other, N);

  for (unsigned int i = 0; i < N; i++) {
    bool d = d_this[i] ^ d_other[i];
    bool e = e_this[i] ^ e_other[i];
    z[i] ^= (x[i] and e) ^ (y[i] and d);
    if (server_num == 0)
      z[i] ^= (d and e);
  }

  delete[] d_this;
  delete[] e_this;
  delete[] d_other;
  delete[] e_other;

  if (do_fork) waitpid(pid, &status, 0);
  return z;
}

// c_{i+1} = c_i xor ((x_i xor c_i) and (y_i xor c_i))
// output z_i = x_i xor y_i xor c_i
// Unused
bool* CorrelatedStore::addBinaryShares(const size_t N,
                                       const size_t* const num_bits,
                                       const bool* const * const x,
                                       const bool* const * const y,
                                       bool* const * const z) {
  bool* carry = new bool[N];
  memset(carry, false, N);
  bool* xi = new bool[N];
  bool* yi = new bool[N];

  size_t max_bits = 0;
  size_t total_bits = 0;
  for (unsigned int i = 0; i < N; i++) {
    max_bits = (num_bits[i] > max_bits ? num_bits[i] : max_bits);
    total_bits += num_bits[i];
  }

  checkBoolTriples(total_bits);

  for (unsigned int j = 0; j < max_bits; j++) {
    size_t idx = 0;
    for (unsigned int i = 0; i < N; i++) {
      if (j >= num_bits[i])
        continue;
      z[i][j] = carry[i] ^ x[i][j] ^ y[i][j];
      xi[idx] = carry[i] ^ x[i][j];
      yi[idx] = carry[i] ^ y[i][j];
      idx++;
    }

    bool* new_carry = multiplyBoolShares(idx, xi, yi);

    idx = 0;
    for (unsigned int i = 0; i < N; i++) {
      if (j >= num_bits[i])
        continue;
      carry[i] ^= new_carry[idx];
      idx++;
    }

    delete[] new_carry;
  }

  delete[] xi;
  delete[] yi;

  return carry;
}

fmpz_t* CorrelatedStore::b2a_daBit_single(const size_t N, const bool* const x) {
  fmpz_t* xp; new_fmpz_array(&xp, N);

  checkDaBits(N);

  bool* v_this = new bool[N];
  for (unsigned int i = 0; i < N; i++) {
    DaBit* dabit = getDaBit();
    v_this[i] = x[i] ^ dabit->b2;

    fmpz_set(xp[i], dabit->bp);
    // consume the daBit
    delete dabit;
  }

  pid_t pid = 0;
  int status = 0;
  if (do_fork) pid = fork();
  if (pid == 0) {
    send_bool_batch(serverfd, v_this, N);

    if (do_fork) exit(EXIT_SUCCESS);
  }
  bool* v_other = new bool[N];
  recv_bool_batch(serverfd, v_other, N);

  for (unsigned int i = 0; i < N; i++) {
    const bool v = v_this[i] ^ v_other[i];

    // [x]_p = v + [b]_p - 2 v [b]_p. Note v only added for one server.
    // So since server_num in {0, 1}, we add it when v = 1
    // Currently, [x]_p is holding [b]_p, which is what we want for v = 0
    if (v) {  // If v = 1, then [x]_p = (0/1) - [b]_p
      fmpz_neg(xp[i], xp[i]);
      fmpz_add_ui(xp[i], xp[i], server_num);
      fmpz_mod(xp[i], xp[i], Int_Modulus);
    }
  }

  delete[] v_this;
  delete[] v_other;

  if (do_fork) waitpid(pid, &status, 0);

  return xp;
}

fmpz_t* CorrelatedStore::b2a_daBit_multi(const size_t N,
                                         const size_t* const num_bits,
                                         const fmpz_t* const x) {
  size_t total_bits = 0;
  for (unsigned int i = 0; i < N; i++)
    total_bits += num_bits[i];

  checkDaBits(total_bits);

  fmpz_t* xp; new_fmpz_array(&xp, N);
  bool* x2 = new bool[total_bits];

  size_t offset = 0;
  for (unsigned int i = 0; i < N; i++) {
    for (unsigned int j = 0; j < num_bits[i]; j++)
      x2[j + offset] = fmpz_tstbit(x[i], j);
    offset += num_bits[i];
  }

  fmpz_t* tmp_xp = b2a_daBit_single(total_bits, x2);

  offset = 0;
  for (unsigned int i = 0; i < N; i++) {
    fmpz_set_ui(xp[i], 0);
    for (unsigned int j = 0; j < num_bits[i]; j++) {
      fmpz_addmul_ui(xp[i], tmp_xp[j + offset], (1 << j));
      fmpz_mod(xp[i], xp[i], Int_Modulus);
    }
    offset += num_bits[i];
  }

  delete[] x2;
  clear_fmpz_array(tmp_xp, total_bits);

  return xp;
}

// Using intsum_ot, multiple bits
fmpz_t* CorrelatedStore::b2a_ot(const size_t num_shares, const size_t num_values,
                                const size_t* const num_bits,
                                const fmpz_t* const x, const size_t mod) {
  uint64_t** x2 = new uint64_t*[num_shares];
  bool* const valid = new bool[num_shares];
  for (unsigned int i = 0; i < num_shares; i++) {
    x2[i] = new uint64_t[num_values];
    valid[i] = true;
    for (unsigned int j = 0; j < num_values; j++) {
      x2[i][j] = fmpz_get_ui(x[i * num_values + j]);
    }
  }

  uint64_t** xp;

  if (server_num == 0) {
    xp = intsum_ot_sender(ot0, x2, valid, num_bits, num_shares, num_values, mod);
  } else {
    xp = intsum_ot_receiver(ot0, x2, num_bits, num_shares, num_values, mod);
  }

  // for consistency, flatten and fmpz_t
  fmpz_t* ans; new_fmpz_array(&ans, num_shares * num_values);

  for (unsigned int i = 0; i < num_shares; i++) {
    for (unsigned int j = 0; j < num_values; j++) {
      fmpz_set_ui(ans[i * num_values + j], xp[i][j]);
    }
    delete[] x2[i];
    delete[] xp[i];
  }
  delete[] x2;
  delete[] valid;
  delete[] xp;

  return ans;
}

// Use b2A via OT on random bit
// Nearly COT, except delta is changing
// random choice and random base, but also random delta matters
DaBit** CorrelatedStore::generateDaBit(const size_t N) {
  DaBit** const dabit = new DaBit*[N];

  emp::PRG prg;
  const size_t mod = fmpz_get_ui(Int_Modulus);

  bool* const b = new bool[N];
  prg.random_bool(b, N);  // random bits

  uint64_t* const x = new uint64_t[N];

  for (unsigned int i = 0; i < N; i++) {
    dabit[i] = new DaBit();
    dabit[i]->b2 = b[i];
  }

  if (server_num == 0) {
    uint64_t* const b0 = new uint64_t[N];
    uint64_t* const b1 = new uint64_t[N];
    for (unsigned int i = 0; i < N; i++) {
      prg.random_data(&b0[i], sizeof(uint64_t));
      b0[i] %= mod;
      b1[i] = (b0[i] + b[i]) % mod;
      x[i] = mod - b0[i];
    }
    ot0->send(b0, b1, N);
    delete[] b0;
    delete[] b1;
  } else {
    ot0->recv(x, b, N);
  }
  for (unsigned int i = 0; i < N; i++) {
    uint64_t bp = (b[i] + 2 * (mod - x[i])) % mod;
    fmpz_set_ui(dabit[i]->bp, bp);
  }

  delete[] b;
  delete[] x;

  return dabit;
}
