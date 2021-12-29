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
  std::queue<const BooleanBeaverTriple* const> new_triples = gen_boolean_beaver_triples(server_num, num_to_make, ot0, ot1);
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
    const DaBit* const * const dabit = generateDaBit(num_to_make);
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

const BooleanBeaverTriple* const CorrelatedStore::getBoolTriple() {
  checkBoolTriples(1);
  const BooleanBeaverTriple* const ans = btriple_store.front();
  btriple_store.pop();
  return ans;
}

const DaBit* const CorrelatedStore::getDaBit() {
  checkDaBits(1);
  const DaBit* const ans = dabit_store.front();
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
    const DaBit* const bit = dabit_store.front();
    dabit_store.pop();
    delete bit;
  }
  while (!btriple_store.empty()) {
    const BooleanBeaverTriple* const triple = btriple_store.front();
    btriple_store.pop();
    delete triple;
  }
}

const bool* const CorrelatedStore::multiplyBoolShares(
    const size_t N, const bool* const x, const bool* const y) {
  bool* const z = new bool[N];

  bool* const d_this = new bool[N];
  bool* const e_this = new bool[N];

  checkBoolTriples(N);
  for (unsigned int i = 0; i < N; i++) {
    const BooleanBeaverTriple* const triple = getBoolTriple();
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
const bool* const CorrelatedStore::addBinaryShares(
    const size_t N, const size_t* const num_bits,
    const bool* const * const x, const bool* const * const y,
    bool* const * const z) {
  bool* const carry = new bool[N];
  memset(carry, false, N);
  bool* const xi = new bool[N];
  bool* const yi = new bool[N];

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

    const bool* const new_carry = multiplyBoolShares(idx, xi, yi);

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

fmpz_t* const CorrelatedStore::b2a_daBit_single(const size_t N, const bool* const x) {
  fmpz_t* xp; new_fmpz_array(&xp, N);

  checkDaBits(N);

  bool* const v_this = new bool[N];
  for (unsigned int i = 0; i < N; i++) {
    const DaBit* const dabit = getDaBit();
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
  bool* const v_other = new bool[N];
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

fmpz_t* const CorrelatedStore::b2a_daBit_multi(
    const size_t N, const size_t* const num_bits, const fmpz_t* const x) {
  size_t total_bits = 0;
  for (unsigned int i = 0; i < N; i++)
    total_bits += num_bits[i];

  checkDaBits(total_bits);

  fmpz_t* xp; new_fmpz_array(&xp, N);
  bool* const x2 = new bool[total_bits];

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
      fmpz_addmul_ui(xp[i], tmp_xp[j + offset], (1ULL << j));
      fmpz_mod(xp[i], xp[i], Int_Modulus);
    }
    offset += num_bits[i];
  }

  delete[] x2;
  clear_fmpz_array(tmp_xp, total_bits);

  return xp;
}

// Using intsum_ot, multiple bits
fmpz_t* const CorrelatedStore::b2a_ot(
    const size_t num_shares, const size_t num_values, const size_t* const num_bits,
    const fmpz_t* const x, const size_t mod) {
  uint64_t** const x2 = new uint64_t*[num_shares];
  bool* const valid = new bool[num_shares];
  for (unsigned int i = 0; i < num_shares; i++) {
    x2[i] = new uint64_t[num_values];
    valid[i] = true;
    for (unsigned int j = 0; j < num_values; j++) {
      x2[i][j] = fmpz_get_ui(x[i * num_values + j]);
    }
  }

  const uint64_t* const * xp;

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

// shares are size 2*N*b = 2 * n
// valid is size N
// buckets are size b
void CorrelatedStore::heavy_ot(
    const size_t N, const size_t b,
    const bool* const shares_x0, const bool* const shares_x1,
    const bool* const valid,
    fmpz_t* const bucket0, fmpz_t* const bucket1) {
  const size_t n = N * b;

  // Step 1: convert x1 as 01 || 11, single bit per
  fmpz_t* shares_p = b2a_daBit_single(2 * n, shares_x1);

  // Step 2: Buckets. Goal is (z + z')(x ^ x')
  // z = servernum - 2y
  // send: (r, r+z) order on x, so (r+xz, r+(1-x)z
  //   and add -r
  // recv: pick using x, get z(x ^ x')

  // 2.1: Setup
  fmpz_t r; fmpz_init(r);
  fmpz_t z; fmpz_init(z);
  fmpz_t tmp; fmpz_init(tmp);
  uint64_t* const data0 = new uint64_t[2*n];
  uint64_t* const data1 = new uint64_t[2*n];
  uint64_t* const received = new uint64_t[2*n];

  for (unsigned int i = 0; i < N; i++) {
    if (!valid[i]) {
      memset(&data0[i*b], 0, b * sizeof(uint64_t));
      memset(&data0[i*b + n], 0, b * sizeof(uint64_t));
      memset(&data1[i*b], 0, b * sizeof(uint64_t));
      memset(&data1[i*b + n], 0, b * sizeof(uint64_t));
      continue;
    }

    for (unsigned int j = 0; j < b; j++) {
      // bucket 0
      size_t idx = i*b + j;
      // rand r, add -r to bucket
      // fmpz_randm(r, seed, Int_Modulus);
      fmpz_zero(r);
      fmpz_sub(bucket0[j], bucket0[j], r);
      fmpz_mod(bucket0[j], bucket0[j], Int_Modulus);
      // z = 1 - 2y
      // std::cout << "b0[" << server_num << ", " << i << ", " << j << "]: subtract r = " << fmpz_get_ui(r) << std::endl;
      // std::cout << "b0[" << server_num << ", " << i << ", " << j << "]: x = " << shares_x0[idx] << std::endl;
      // std::cout << "b0[" << server_num << ", " << i << ", " << j << "]: y = " << fmpz_get_ui(shares_p[idx]) << std::endl;
      fmpz_set_si(z, server_num);
      fmpz_submul_si(z, shares_p[idx], 2);
      fmpz_mod(z, z, Int_Modulus);
      // std::cout << "b0[" << server_num << ", " << i << ", " << j << "]: z = " << fmpz_get_ui(z) << std::endl;
      // 0 = r + xz
      fmpz_set(tmp, r); fmpz_addmul_ui(tmp, z, shares_x0[idx]);
      fmpz_mod(tmp, tmp, Int_Modulus);
      data0[idx] = fmpz_get_ui(tmp);
      // 1 = r + (1-x)z
      fmpz_set(tmp, r); fmpz_addmul_ui(tmp, z, 1 - shares_x0[idx]);
      fmpz_mod(tmp, tmp, Int_Modulus);
      data1[idx] = fmpz_get_ui(tmp);
      // std::cout << "b0[" << server_num << ", " << i << ", " << j << "]: send (" << data0[idx] << ", " << data1[idx] << ")\n";

      // bucket 1, idx += n
      idx = i*b + j + n;
      // rand r, add -r to bucket
      // fmpz_randm(r, seed, Int_Modulus);
      fmpz_zero(r);
      fmpz_sub(bucket1[j], bucket1[j], r);
      fmpz_mod(bucket1[j], bucket1[j], Int_Modulus);
      // z = 1 - 2y
      // std::cout << "b1[" << server_num << ", " << i << ", " << j << "]: subtract r = " << fmpz_get_ui(r) << std::endl;
      // std::cout << "b1[" << server_num << ", " << i << ", " << j << "]: x = " << shares_x0[idx] << std::endl;
      // std::cout << "b1[" << server_num << ", " << i << ", " << j << "]: y = " << fmpz_get_ui(shares_p[idx]) << std::endl;
      fmpz_set_si(z, server_num);
      fmpz_submul_si(z, shares_p[idx], 2);
      fmpz_mod(z, z, Int_Modulus);
      // std::cout << "b1[" << server_num << ", " << i << ", " << j << "]: z = " << fmpz_get_ui(z) << std::endl;
      // 0 = r + xz
      fmpz_set(tmp, r); fmpz_addmul_ui(tmp, z, shares_x0[idx]);
      fmpz_mod(tmp, tmp, Int_Modulus);
      data0[idx] = fmpz_get_ui(tmp);
      // 1 = r + (1-x)z
      fmpz_set(tmp, r); fmpz_addmul_ui(tmp, z, 1 - shares_x0[idx]);
      fmpz_mod(tmp, tmp, Int_Modulus);
      data1[idx] = fmpz_get_ui(tmp);
      // std::cout << "b1[" << server_num << ", " << i << ", " << j << "]: send (" << data0[idx] << ", " << data1[idx] << ")\n";
    }
  }

  // OT swap
  pid_t pid = 0;
  int status = 0;
  // OT forking currently seems bugged. Disable for now. 
  const bool do_ot_fork = false;
  if (do_ot_fork) {
    pid = fork();
    if (pid == 0) {
      (server_num == 0 ? ot0 : ot1)->send(data0, data1, 2 * n);
      exit(EXIT_SUCCESS);
    }
    (server_num == 0 ? ot1 : ot0)->recv(received, shares_x0, 2 * n);
  } else {
    if (server_num == 0) {
      ot0->send(data0, data1, 2 * n);
      ot1->recv(received, shares_x0, 2 * n);
    } else {
      ot0->recv(received, shares_x0, 2 * n);
      ot1->send(data0, data1, 2 * n);
    }
  }

  // add received to buckets
  for (unsigned int i = 0; i < N; i++) {
    for (unsigned int j = 0; j < b; j++) {
      size_t idx = i * b + j;
      // std::cout << "(" << i << ", " << j << ")" << std::endl;
      fmpz_add_ui(bucket0[j], bucket0[j], received[idx]);
      fmpz_mod(bucket0[j], bucket0[j], Int_Modulus);
      // std::cout << " b0[" << server_num << ", " << i << ", " << j << "]: add recv[" << idx << "] = " << received[idx] << std::endl;

      idx = (i * b + j) + n;
      fmpz_add_ui(bucket1[j], bucket1[j], received[idx]);
      fmpz_mod(bucket1[j], bucket1[j], Int_Modulus);
      // std::cout << " b1[" << server_num << ", " << i << ", " << j << "]: add recv[" << idx << "] = " << received[idx] << std::endl;
    }
  }

  fmpz_clear(r);
  fmpz_clear(z);
  fmpz_clear(tmp);
  delete[] data0;
  delete[] data1;
  delete[] received;

  if (do_ot_fork) waitpid(pid, &status, 0);
}

// Use b2A via OT on random bit
// Nearly COT, except delta is changing
// random choice and random base, but also random delta matters
const DaBit* const * const CorrelatedStore::generateDaBit(const size_t N) {
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

// TODO: do better. Currently just in clear.
// [x] > c for known c, shares [x]
// Treats > N/2 as negative
// Maybe also return int (+/-/0) for equality?
bool* CorrelatedStore::cmp_c(const size_t N,
                             const fmpz_t* const x,
                             const fmpz_t* const c) {
  bool* ans = new bool[N];

  fmpz_t half; fmpz_init(half); fmpz_cdiv_q_ui(half, Int_Modulus, 2);

  // For now, just do in clear. TODO: do securely.
  if (server_num == 0) {
    fmpz_t* x2; new_fmpz_array(&x2, N);
    fmpz_t* c2; new_fmpz_array(&c2, N);
    recv_fmpz_batch(serverfd, x2, N);
    for (unsigned int i = 0; i < N; i++) {
      fmpz_add(x2[i], x2[i], x[i]);
      fmpz_mod(x2[i], x2[i], Int_Modulus);
      if (fmpz_cmp(x2[i], half) > 0) {  // > N/2, so negative
        fmpz_sub(x2[i], x2[i], Int_Modulus);
      }
      if (fmpz_cmp(c[i], half) > 0) {
        fmpz_sub(c2[i], c[i], Int_Modulus);
      } else {
        fmpz_set(c2[i], c[i]);
      }
      ans[i] = fmpz_cmp(x2[i], c2[i]) < 0;
    }
    send_bool_batch(serverfd, ans, N);
  } else {
    send_fmpz_batch(serverfd, x, N);
    recv_bool_batch(serverfd, ans, N);
  }

  fmpz_clear(half);
  return ans;
}

// Just x < y => (x - y) < 0
bool* CorrelatedStore::cmp(const size_t N,
                           const fmpz_t* const x,
                           const fmpz_t* const y) {
  fmpz_t* diff; new_fmpz_array(&diff, N);
  fmpz_t* zeros; new_fmpz_array(&zeros, N);  // zeroed out by default

  fmpz_t half; fmpz_init(half); fmpz_cdiv_q_ui(half, Int_Modulus, 2);
  fmpz_t tmp; fmpz_init(tmp);

  for (unsigned int i = 0; i < N; i++) {
    fmpz_sub(diff[i], x[i], y[i]);
    fmpz_mod(diff[i], diff[i], Int_Modulus);
    // std::cout << "diff_" << server_num << "[" << i << "] = " << fmpz_get_ui(diff[i]) << std::endl;
  }

  bool* ans = cmp_c(N, diff, zeros);
  fmpz_clear(tmp);
  fmpz_clear(half);
  clear_fmpz_array(diff, N);
  clear_fmpz_array(zeros, N);
  return ans;
}

void CorrelatedStore::abs(const size_t N, const fmpz_t* const x, fmpz_t* out) {
  fmpz_t* zeros; new_fmpz_array(&zeros, N);
  bool* sign = cmp_c(N, x, zeros);
  for (unsigned int i = 0; i < N; i++) {
    if (sign[i]) {  // negative, flip
      fmpz_sub(out[i], Int_Modulus, x[i]);
    } else {
      fmpz_set(out[i], x[i]);
    }
  }
}

bool* CorrelatedStore::abs_cmp(const size_t N,
                               const fmpz_t* const x, const fmpz_t* const y) {
  fmpz_t* merge; new_fmpz_array(&merge, 2*N);
  fmpz_t* merge_abs; new_fmpz_array(&merge_abs, 2*N);
  for (unsigned int i = 0; i < N; i++) {
    fmpz_set(merge[i], x[i]);
    fmpz_set(merge[i+N], y[i]);
  }
  abs(2*N, merge, merge_abs);
  clear_fmpz_array(merge, 2*N);

  fmpz_t* x2; new_fmpz_array(&x2, N);
  fmpz_t* y2; new_fmpz_array(&y2, N);
  for (unsigned int i = 0; i < N; i++) {
    fmpz_set(x2[i], merge_abs[i]);
    fmpz_set(y2[i], merge_abs[i+N]);
  }
  clear_fmpz_array(merge_abs, N);
  bool* ans = cmp(N, x2, y2);
  clear_fmpz_array(x2, N);
  clear_fmpz_array(y2, N);
  return ans;
}
