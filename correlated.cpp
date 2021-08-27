#include "correlated.h"

#include <sys/wait.h>

#include <iostream>

#include "constants.h"
#include "net_share.h"
#include "ot.h"
#include "utils.h"

void CorrelatedStore::addBoolTriples(const size_t n) {
  auto start = clock_start();
  const size_t num_to_make = (n > bool_batch_size ? n : bool_batch_size);
  std::cout << "adding booltriples: " << num_to_make << std::endl;
  std::queue<BooleanBeaverTriple*> new_triples = gen_boolean_beaver_triples(server_num, num_to_make, ot0, ot1);
  for (unsigned int i = 0; i < num_to_make; i++) {
    btriple_store.push(new_triples.front());
    new_triples.pop();
  }
  std::cout << "addBoolTriples timing : " << sec_from(start) << std::endl;
}

void CorrelatedStore::addTriples(const size_t n) {
  auto start = clock_start();
  const size_t num_to_make = (n > batch_size ? n : batch_size);
  std::cout << "adding triples: " << num_to_make << std::endl;
  if (triple_gen) {  // not null pointer
    std::vector<BeaverTriple*> new_triples = triple_gen->generateTriples(num_to_make);
    for (unsigned int i = 0; i < num_to_make; i++)
      atriple_store.push(new_triples[i]);
  } else {
    std::cout << "Using lazy beaver triples" << std::endl;
    // std::cout << "Using OT beaver triples" << std::endl;
    for (unsigned int i = 0; i < num_to_make; i++) {
      BeaverTriple* triple = generate_beaver_triple_lazy(serverfd, server_num);
      // BeaverTriple* triple = generate_beaver_triple(serverfd, server_num, ot0, ot1);
      atriple_store.push(triple);
    }
  }
  std::cout << "addTriples timing : " << sec_from(start) << std::endl;
}

void CorrelatedStore::addDaBits(const size_t n) {
  auto start = clock_start();
  const size_t num_to_make = (n > batch_size ? n : batch_size);
  std::cout << "adding dabits: " << num_to_make << std::endl;
  if (!lazy) {
    DaBit** dabit = generateDaBit(num_to_make);
    for (unsigned int i = 0; i < num_to_make; i++)
      dabit_store.push(dabit[i]);
    delete[] dabit;
  } else {
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
      delete[] dabit;
    }
    else {
      recv_DaBit_batch(serverfd, dabit, num_to_make);
      for (unsigned int i = 0; i < num_to_make; i++)
        dabit_store.push(dabit[i]);
      delete[] dabit;
    }
  }
  std::cout << "addDaBits timing : " << sec_from(start) << std::endl;
}

void CorrelatedStore::addEdaBits(const size_t num_bits, const size_t n) {
  auto start = clock_start();
  const size_t num_to_make = (n > batch_size ? n : batch_size);
  std::cout << "adding " << num_bits << " bit edabits: " << num_to_make << std::endl;
  if (lazy) {
    if (server_num == 0) {  // make on server 0
      std::cout << "Making lazy edabits" << std::endl;
      EdaBit* other_edabit = new EdaBit(num_bits);
      for (unsigned int i = 0; i < num_to_make; i++) {
        EdaBit* edabit = new EdaBit(num_bits);
        makeLocalEdaBit(edabit, other_edabit, num_bits);
        send_EdaBit(serverfd, other_edabit, num_bits);  // TODO: batch
        if (num_bits == nbits)
          edabit_store.push(edabit);
        else if (num_bits == 2 * nbits)
          edabit_store_2.push(edabit);
      }
      delete other_edabit;
    } else {
      for (unsigned int i = 0; i < num_to_make; i++) {
        EdaBit* edabit = new EdaBit(num_bits);
        recv_EdaBit(serverfd, edabit, num_bits);
        if (num_bits == nbits)
          edabit_store.push(edabit);
        else if (num_bits == 2 * nbits)
          edabit_store_2.push(edabit);
      }
    }
    std::cout << "lazy addEdaBits timing : " << sec_from(start) << std::endl;
    return;
  }

  EdaBit** edabit = generateEdaBit(num_to_make, num_bits);
  for (unsigned int i = 0; i < num_to_make; i++) {
    if (num_bits == nbits)
      edabit_store.push(edabit[i]);
    else if (num_bits == 2 * nbits)
      edabit_store_2.push(edabit[i]);
  }
  delete[] edabit;

  std::cout << "addEdaBits timing : " << sec_from(start) << std::endl;
}

void CorrelatedStore::checkBoolTriples(const size_t n) { 
  if (btriple_store.size() < n) addBoolTriples(n - btriple_store.size());
}

void CorrelatedStore::checkTriples(const size_t n, const bool always) { 
  if ((!lazy or always) and atriple_store.size() < n) addTriples(n - atriple_store.size());
}

void CorrelatedStore::checkDaBits(const size_t n) { 
  if (dabit_store.size() < n) addDaBits(n - dabit_store.size());
}

void CorrelatedStore::checkEdaBits(const size_t num_bits, const size_t n) { 
  if (num_bits == nbits) {
    if (edabit_store.size() < n) addEdaBits(num_bits, n - edabit_store.size());
  } else if (num_bits == 2 * nbits) {
    if (edabit_store_2.size() < n) addEdaBits(num_bits, n - edabit_store_2.size());
  } else {
    std::cerr << "Don't support " << num_bits << " bit edabits" << std::endl;
    std::cerr << "Only " << nbits << " or " << 2 * nbits << " is supported" << std::endl;
    exit(EXIT_FAILURE);
  }
}

BooleanBeaverTriple* CorrelatedStore::getBoolTriple() {
  checkBoolTriples(1);
  BooleanBeaverTriple* ans = btriple_store.front();
  btriple_store.pop();
  return ans;
}

BeaverTriple* CorrelatedStore::getTriple() {
  checkTriples(1);
  BeaverTriple* ans = atriple_store.front();
  atriple_store.pop();
  return ans;
}

DaBit* CorrelatedStore::getDaBit() {
  checkDaBits(1);
  DaBit* ans = dabit_store.front();
  dabit_store.pop();
  return ans;
}

EdaBit* CorrelatedStore::getEdaBit(const size_t num_bits) {
  checkEdaBits(num_bits, 1);
  EdaBit* ans;
  if (num_bits == nbits) {
    ans = edabit_store.front();
    edabit_store.pop();
    return ans;
  } else if (num_bits == 2 * nbits) {
    ans = edabit_store_2.front();
    edabit_store_2.pop();
  } else {
    std::cerr << "Cannot get a " << num_bits << " bit edabits" << std::endl;
    std::cerr << "Only " << nbits << " or " << 2 * nbits << " is supported" << std::endl;
    exit(EXIT_FAILURE);
  }

  return ans;
}

void CorrelatedStore::printSizes() {
  std::cout << "Current store sizes:" << std::endl;
  std::cout << "       EdaBits: " << edabit_store.size() << std::endl;
  std::cout << "     EdaBits 2: " << edabit_store_2.size() << std::endl;
  std::cout << "        Dabits: " << dabit_store.size() << std::endl;
  std::cout << " Arith Triples: " << atriple_store.size() << std::endl;
  std::cout << " Bool  Triples: " << btriple_store.size() << std::endl;
}

void CorrelatedStore::maybeUpdate(const bool using_eda) {
  std::cout << "precomputing using " << (using_eda?"e":"") << "dabits..." << std::endl;
  auto start = clock_start();

  // If making extra
  const size_t extra = over_precompute ? 1 : 0;
  // Make top level if stores not enough
  const bool make_eda = using_eda && (edabit_store.size() < batch_size / 2);
  const bool make_eda2 = using_eda && (edabit_store_2.size() < batch_size / 2);
  const bool make_da = !using_eda && (dabit_store.size() < batch_size * nbits / 2);
  // Determine how much of each to make
  const size_t da_target = (using_eda 
                            ? (make_eda + make_eda2 + extra) 
                            : (2 * make_da * nbits));
  const size_t atrip_target = da_target + extra;
  const size_t btrip_target = using_eda ? (2 * (make_eda + 2 * make_eda2) + extra) : 0;

  if (btriple_store.size() < btrip_target * bool_batch_size)
      addBoolTriples(btrip_target * bool_batch_size);

  if (!lazy) {
    if (atriple_store.size() < atrip_target * batch_size)
      addTriples(atrip_target * batch_size);
  }
  if (!lazy or !using_eda) {
    if (dabit_store.size() < da_target * batch_size)
      addDaBits(da_target * batch_size);
  }

  if (make_eda) addEdaBits(nbits);
  if (make_eda2) addEdaBits(2 * nbits);

  printSizes();

  std::cout << "precompute timing : " << sec_from(start) << std::endl;
}

CorrelatedStore::~CorrelatedStore() {
  while (!edabit_store.empty()) {
    EdaBit* bit = edabit_store.front();
    edabit_store.pop();
    delete bit;
  }
  while (!edabit_store_2.empty()) {
    EdaBit* bit = edabit_store_2.front();
    edabit_store_2.pop();
    delete bit;
  }
  while (!dabit_store.empty()) {
    DaBit* bit = dabit_store.front();
    dabit_store.pop();
    delete bit;
  }
  while (!atriple_store.empty()) {
    BeaverTriple* triple = atriple_store.front();
    atriple_store.pop();
    delete triple;
  }
  while (!btriple_store.empty()) {
    BooleanBeaverTriple* triple = btriple_store.front();
    btriple_store.pop();
    delete triple;
  }
  if (triple_gen)
    delete triple_gen;
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

fmpz_t* CorrelatedStore::multiplyArithmeticShares(const size_t N,
                                                  const fmpz_t* const x,
                                                  const fmpz_t* const y) {
  fmpz_t* z; new_fmpz_array(&z, N);

  fmpz_t* d; new_fmpz_array(&d, N);
  fmpz_t* e; new_fmpz_array(&e, N);
  fmpz_t* d_other; new_fmpz_array(&d_other, N);
  fmpz_t* e_other; new_fmpz_array(&e_other, N);

  checkTriples(N, true);

  for (unsigned int i = 0; i < N; i++) {
    BeaverTriple* triple = getTriple();

    fmpz_sub(d[i], x[i], triple->A);  // [d] = [x] - [a]
    fmpz_mod(d[i], d[i], Int_Modulus);
    fmpz_sub(e[i], y[i], triple->B);  // [e] = [y] - [b]
    fmpz_mod(e[i], e[i], Int_Modulus);

    fmpz_set(z[i], triple->C);

    // consume the triple
    delete triple;
  }

  // Spawn a child to do the sending, so that can recieve at the same time
  pid_t pid = 0;
  int status = 0;
  if (do_fork) pid = fork();
  if (pid == 0) {
    send_fmpz_batch(serverfd, d, N);
    send_fmpz_batch(serverfd, e, N);
    if (do_fork) exit(EXIT_SUCCESS);
  }
  recv_fmpz_batch(serverfd, d_other, N);
  recv_fmpz_batch(serverfd, e_other, N);

  for (unsigned int i = 0; i < N; i++) {
    fmpz_add(d[i], d[i], d_other[i]);  // x - a
    fmpz_mod(d[i], d[i], Int_Modulus);
    fmpz_add(e[i], e[i], e_other[i]);  // y - b
    fmpz_mod(e[i], e[i], Int_Modulus);

    // [xy] = [c] + [x] e + [y] d - de
    // Is it more efficient using a tmp for product, and modding more?
    fmpz_addmul(z[i], x[i], e[i]);
    fmpz_addmul(z[i], y[i], d[i]);
    if (server_num == 0)
      fmpz_submul(z[i], d[i], e[i]);
    fmpz_mod(z[i], z[i], Int_Modulus);
  }

  clear_fmpz_array(d_other, N);
  clear_fmpz_array(e_other, N);
  clear_fmpz_array(d, N);
  clear_fmpz_array(e, N);

  // Wait for send child to finish
  if (do_fork) waitpid(pid, &status, 0);

  return z;
}

// c_{i+1} = c_i xor ((x_i xor c_i) and (y_i xor c_i))
// output z_i = x_i xor y_i xor c_i
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

fmpz_t* CorrelatedStore::b2a_edaBit(const size_t N,
                                    const size_t* const num_bits,
                                    const fmpz_t* const x) {
  size_t num_n = 0, num_2n = 0;
  for (unsigned int i = 0; i < N; i++) {
    if (num_bits[i] == nbits) num_n += 1;
    if (num_bits[i] == 2 * nbits) num_2n += 1;
  }
  const size_t eda_to_make = (edabit_store.size() < num_n) ? num_n - edabit_store.size() : 0;
  const size_t eda2_to_make = (edabit_store_2.size() < num_2n) ? num_2n - edabit_store_2.size() : 0;
  const size_t da_target = eda_to_make + eda2_to_make;
  const size_t bool_target = (1 + !lazy) * nbits * (eda_to_make + 2 * eda2_to_make);

  checkBoolTriples(bool_target);
  checkTriples(da_target);
  checkDaBits(da_target);
  checkEdaBits(nbits, num_n);
  checkEdaBits(2 * nbits, num_2n);

  fmpz_t* xp; new_fmpz_array(&xp, N);

  bool** x2 = new bool*[N];
  bool** b = new bool*[N];
  bool** xr = new bool*[N];
  fmpz_t* ebit_r; new_fmpz_array(&ebit_r, N);

  for (unsigned int i = 0; i < N; i++) {
    x2[i] = new bool[num_bits[i]];
    b[i] = new bool[num_bits[i]];
    xr[i] = new bool[num_bits[i] + 1];
    EdaBit* edabit = getEdaBit(num_bits[i]);
    for (unsigned int j = 0; j < num_bits[i]; j++) {
      // Convert x2 to bool array
      x2[i][j] = fmpz_tstbit(x[i], j);
      b[i][j] = edabit->b[j];
    }
    fmpz_set(ebit_r[i], edabit->r);

    // consume edabit
    delete edabit;
  }

  // [x + r]_2 = [x]_2 + [r]_2 via circuit
  bool* carry = addBinaryShares(N, num_bits, x2, b, xr);

  for (unsigned int i = 0; i < N; i++) {
    xr[i][num_bits[i]] = carry[i];

    fmpz_from_bool_array(xp[i], xr[i], num_bits[i] + 1);  // [x + r]_2

    delete[] x2[i];
    delete[] b[i];
    delete[] xr[i];
  }
  delete[] carry;
  delete[] x2;
  delete[] b;
  delete[] xr;

  // reveal x + r, convert to mod p shares
  if (server_num == 0) {
    fmpz_t* xr_other; new_fmpz_array(&xr_other, N);
    recv_fmpz_batch(serverfd, xr_other, N);  // get other [x + r]_2
    for (unsigned int i = 0; i < N; i++) {
      fmpz_xor(xp[i], xp[i], xr_other[i]);  // real x + r
      fmpz_randm(xr_other[i], seed, Int_Modulus);  // make other [x + r]_p
      fmpz_sub(xp[i], xp[i], xr_other[i]);
      fmpz_mod(xp[i], xp[i], Int_Modulus);  // This [x + r]_p
    }
    send_fmpz_batch(serverfd, xr_other, N);
    clear_fmpz_array(xr_other, N);
  } else {
    send_fmpz_batch(serverfd, xp, N);
    recv_fmpz_batch(serverfd, xp, N);
  }

  // [x]_p = [x+r]_p - [r]_p
  for (unsigned int i = 0; i < N; i++) {
    fmpz_sub(xp[i], xp[i], ebit_r[i]);
    fmpz_mod(xp[i], xp[i], Int_Modulus);
  }

  clear_fmpz_array(ebit_r, N);

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
  }
  delete[] x2;
  delete[] valid;

  return ans;
}

DaBit** CorrelatedStore::generateDaBit(const size_t N) {
  DaBit** const dabit = new DaBit*[N];

  DaBit** const dabit0 = new DaBit*[N];
  DaBit** const dabit1 = new DaBit*[N];
  fmpz_t* x; new_fmpz_array(&x, N);
  fmpz_t* y; new_fmpz_array(&y, N);

  for (unsigned int i = 0; i < N; i++) {
    dabit[i] = new DaBit();

    // Create local
    dabit0[i] = new DaBit();
    dabit1[i] = new DaBit();
    makeLocalDaBit(dabit0[i], dabit1[i]);
  }
  // Exchange
  pid_t pid = 0;
  int status = 0;
  if (do_fork) pid = fork();
  if (pid == 0) {
    send_DaBit_batch(serverfd, server_num == 0 ? dabit1 : dabit0, N);
    if (do_fork) exit(EXIT_SUCCESS);
  }

  recv_DaBit_batch(serverfd, server_num == 0 ? dabit1 : dabit0, N);
  for (unsigned int i = 0; i < N; i++) {
    // Xor boolean shares
    dabit[i]->b2 = dabit0[i]->b2 ^ dabit1[i]->b2;

    fmpz_set(x[i], dabit0[i]->bp);
    fmpz_set(y[i], dabit1[i]->bp);

    delete dabit0[i];
    delete dabit1[i];
  }
  delete[] dabit0;
  delete[] dabit1;

  if (do_fork) waitpid(pid, &status, 0);

  // Xor Arithmetic shares, using a xor b = a + b - 2ab
  fmpz_t* z = multiplyArithmeticShares(N, x, y);

  for (unsigned int i = 0; i < N; i++) {
    fmpz_add(dabit[i]->bp, x[i], y[i]);
    fmpz_submul_ui(dabit[i]->bp, z[i], 2);
    fmpz_mod(dabit[i]->bp, dabit[i]->bp, Int_Modulus);
  }

  clear_fmpz_array(x, N);
  clear_fmpz_array(y, N);
  clear_fmpz_array(z, N);

  return dabit;
}

EdaBit** CorrelatedStore::generateEdaBit(const size_t N, const size_t num_bits) {
  EdaBit** edabit = new EdaBit*[N];

  EdaBit** edabit0 = new EdaBit*[N];
  EdaBit** edabit1 = new EdaBit*[N];
  bool** b = new bool*[N];
  bool** b0 = new bool*[N];
  bool** b1 = new bool*[N];

  for (unsigned int i = 0; i < N; i++) {
    edabit[i] = new EdaBit(num_bits);

    // Create local
    edabit0[i] = new EdaBit(num_bits);
    edabit1[i] = new EdaBit(num_bits);
    makeLocalEdaBit(edabit0[i], edabit1[i], num_bits);
  }
  // Exchange
  pid_t pid = 0;
  int status = 0;
  if (do_fork) pid = fork();
  if (pid == 0) {
    send_EdaBit_batch(serverfd, server_num == 0 ? edabit1 : edabit0, num_bits, N);
    if (do_fork) exit(EXIT_SUCCESS);
  }

  recv_EdaBit_batch(serverfd, server_num == 0 ? edabit1 : edabit0, num_bits, N);
  for (unsigned int i = 0; i < N; i++) {
    // Add arithmetic shares
    fmpz_add(edabit[i]->r, edabit0[i]->r, edabit1[i]->r);
    fmpz_mod(edabit[i]->r, edabit[i]->r, Int_Modulus);

    b[i] = new bool[num_bits];
    b0[i] = new bool[num_bits];
    b1[i] = new bool[num_bits];
    memcpy(b0[i], edabit0[i]->b, num_bits * sizeof(bool));
    memcpy(b1[i], edabit1[i]->b, num_bits * sizeof(bool));

    delete edabit0[i];
    delete edabit1[i];
  }
  delete[] edabit0;
  delete[] edabit1;

  if (do_fork) waitpid(pid, &status, 0);

  // Add binary shares via circuit
  size_t* bits_arr = new size_t[N];
  for (unsigned int i = 0; i < N; i++)
    bits_arr[i] = num_bits;
  bool* carry = addBinaryShares(N, bits_arr, b0, b1, b);
  delete[] bits_arr;

  fmpz_t tmp; fmpz_init(tmp);
  fmpz_from_bool_array(tmp, b[0], num_bits);

  // Convert carry to arithmetic [carry]_p
  fmpz_t* carry_p = b2a_daBit_single(N, carry);
  delete[] carry;

  // Subtract out 2^n * [carry]_p from r
  fmpz_t pow; fmpz_init_set_ui(pow, 2);
  fmpz_pow_ui(pow, pow, num_bits);
  for (unsigned int i = 0; i < N; i++) {
    fmpz_submul(edabit[i]->r, carry_p[i], pow);
    fmpz_mod(edabit[i]->r, edabit[i]->r, Int_Modulus);

    memcpy(edabit[i]->b, b[i], num_bits * sizeof(bool));

    delete[] b[i];
    delete[] b0[i];
    delete[] b1[i];
  }
  delete[] b;
  delete[] b0;
  delete[] b1;
  clear_fmpz_array(carry_p, N);
  fmpz_clear(pow);

  return edabit;
}
