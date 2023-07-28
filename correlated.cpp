#include "correlated.h"

#include <sys/wait.h>

#include <iostream>

#include "constants.h"
#include "net_share.h"
#include "ot.h"
#include "utils.h"

const DaBit* const CorrelatedStore::get_DaBit() {
  check_DaBits(1);
  const DaBit* const ans = dabit_store.front();
  dabit_store.pop();
  return ans;
}

const BeaverTriple* const CorrelatedStore::get_Triple() {
  check_Triples(1);
  const BeaverTriple* const ans = atriple_store.front();
  atriple_store.pop();
  return ans;
}

const BooleanBeaverTriple* const CorrelatedStore::get_BoolTriple() {
  check_BoolTriples(1);
  const BooleanBeaverTriple* const ans = btriple_store.front();
  btriple_store.pop();
  return ans;
}

void CorrelatedStore::b2a_single_setup(
    const size_t N, const bool* const x, fmpz_t* const xp, bool* const v) {
  for (unsigned int i = 0; i < N; i++) {
    const DaBit* const dabit = get_DaBit();
    v[i] = x[i] ^ dabit->b2;

    fmpz_set(xp[i], dabit->bp);
    // consume the daBit
    delete dabit;
  }
}

void CorrelatedStore::b2a_single_finish(
    const size_t N, fmpz_t* const xp,
    const bool* const v, const bool* const v_other) {
  for (unsigned int i = 0; i < N; i++) {
    // [x]_p = v + [b]_p - 2 v [b]_p. Note v only added for one server.
    // So since server_num in {0, 1}, we add it when v = 1
    // Currently, [x]_p is holding [b]_p, which is what we want for v = 0
    if (v[i] ^ v_other[i]) {  // If v = 1, then [x]_p = (0/1) - [b]_p
      fmpz_mod_neg(xp[i], xp[i], mod_ctx);
      fmpz_mod_add_ui(xp[i], xp[i], server_num, mod_ctx);
    }
  }
}

int CorrelatedStore::b2a_single(const size_t N, const bool* const x, fmpz_t* const xp) {
  int sent_bytes = 0;

  check_DaBits(N);

  bool* const v = new bool[N];
  bool* const v_other = new bool[N];

  b2a_single_setup(N, x, xp, v);

  sent_bytes += send_bool_batch(serverfd, v, N);
  recv_bool_batch(serverfd, v_other, N);

  b2a_single_finish(N, xp, v, v_other);

  delete[] v;
  delete[] v_other;

  return sent_bytes;
}

void CorrelatedStore::b2a_multi_setup(
    const size_t N, const size_t total_bits, const size_t* const num_bits,
    const fmpz_t* const x, fmpz_t* const flat_xp, bool* const v) {
  bool* const x2 = new bool[total_bits];

  size_t offset = 0;
  for (unsigned int i = 0; i < N; i++) {
    for (unsigned int j = 0; j < num_bits[i]; j++)
      x2[j + offset] = fmpz_tstbit(x[i], j);
    offset += num_bits[i];
  }

  b2a_single_setup(total_bits, x2, flat_xp, v);
  delete[] x2;
}

void CorrelatedStore::b2a_multi_finish(
    const size_t N, const size_t total_bits, const size_t* const num_bits,
    fmpz_t* const xp, fmpz_t* const flat_xp,
    const bool* const v, const bool* const v_other
  ) {

  b2a_single_finish(total_bits, flat_xp, v, v_other);

  size_t offset = 0;
  for (unsigned int i = 0; i < N; i++) {
    fmpz_set_ui(xp[i], 0);
    for (unsigned int j = 0; j < num_bits[i]; j++) {
      fmpz_mod_addmul_ui(xp[i], flat_xp[j + offset], (1ULL << j), mod_ctx);
    }
    offset += num_bits[i];
  }
}

int CorrelatedStore::b2a_multi(
    const size_t N, const size_t* const num_bits,
    const fmpz_t* const x, fmpz_t* const xp) {
  int sent_bytes = 0;

  size_t total_bits = 0;
  for (unsigned int i = 0; i < N; i++)
    total_bits += num_bits[i];

  check_DaBits(total_bits);

  fmpz_t* flat_xp; new_fmpz_array(&flat_xp, total_bits);
  bool* const v = new bool[total_bits];
  bool* const v_other = new bool[total_bits];

  b2a_multi_setup(N, total_bits, num_bits, x, flat_xp, v);

  sent_bytes += send_bool_batch(serverfd, v, total_bits);
  recv_bool_batch(serverfd, v_other, total_bits);

  b2a_multi_finish(N, total_bits, num_bits, xp, flat_xp, v, v_other);

  delete[] v;
  delete[] v_other;
  clear_fmpz_array(flat_xp, total_bits);

  return sent_bytes;
}

void CorrelatedStore::multiply_BoolShares_setup(
    const size_t N, const bool* const x, const bool* const y, bool* const z,
    bool* const de) {
  for (unsigned int i = 0; i < N; i++) {
    const BooleanBeaverTriple* const triple = get_BoolTriple();
    de[i] = x[i] ^ triple->a;
    de[i + N] = y[i] ^ triple->b;
    z[i] = triple->c;
    delete triple;
  }
}

void CorrelatedStore::multiply_BoolShares_finish(
    const size_t N, const bool* const x, const bool* const y, bool* const z,
    const bool* const de, const bool* const de_other) {

    for (unsigned int i = 0; i < N; i++) {
      bool d = de[i] ^ de_other[i];
      bool e = de[i + N] ^ de_other[i + N];
      z[i] ^= (x[i] and e) ^ (y[i] and d);
      if (server_num == 0)
        z[i] ^= (d and e);
    }
}

// 1 bool triple per bit
int CorrelatedStore::multiply_BoolShares(
    const size_t N, const bool* const x, const bool* const y, bool* const z) {
  int sent_bytes = 0;

  check_BoolTriples(2 * N);

  bool* de = new bool[2 * N];
  bool* de_other = new bool[2 * N];

  multiply_BoolShares_setup(N, x, y, z, de);

  sent_bytes += send_bool_batch(serverfd, de, 2 * N);
  recv_bool_batch(serverfd, de_other, 2 * N);

  multiply_BoolShares_finish(N, x, y, z, de, de_other);

  delete[] de;

  return sent_bytes;
}

// 1 triple per bit
int CorrelatedStore::multiply_ArithmeticShares(
    const size_t N, const fmpz_t* const x, const fmpz_t* const y,
    fmpz_t* const z) {
  int sent_bytes = 0;

  fmpz_t* de; new_fmpz_array(&de, 2 * N);

  check_Triples(N);

  for (unsigned int i = 0; i < N; i++) {
    const BeaverTriple* const triple = get_Triple();

    fmpz_mod_sub(de[i], x[i], triple->A, mod_ctx);  // [d] = [x] - [a]
    fmpz_mod_sub(de[i+N], y[i], triple->B, mod_ctx);  // [e] = [y] - [b]

    fmpz_set(z[i], triple->C);

    // consume the triple
    delete triple;
  }

  sent_bytes += reveal_fmpz_batch(serverfd, de, 2 * N);

  for (unsigned int i = 0; i < N; i++) {
    // [xy] = [c] + [x] e + [y] d - de
    // Is it more efficient using a tmp for product, and modding more?
    fmpz_mod_addmul(z[i], x[i], de[i+N], mod_ctx);
    fmpz_mod_addmul(z[i], y[i], de[i], mod_ctx);
    if (server_num == 0)
      fmpz_mod_submul(z[i], de[i], de[i+N], mod_ctx);
  }

  clear_fmpz_array(de, 2 * N);

  return sent_bytes;
}

void CorrelatedStore::multiply_BoolShares_cross_setup(
    const size_t N, const size_t a, const size_t b,
    const bool* x, const bool* y, bool* const z,
    bool* const x_ext, bool* const y_ext, bool* const de) {

  for (unsigned int i = 0; i < N * a; i++) {
    memset(&x_ext[i * b], x[i], b);
    memcpy(&y_ext[i * b], &y[(i / a) * b], b);
  }

  multiply_BoolShares_setup(N * a * b, x_ext, y_ext, z, de);
}

int CorrelatedStore::multiply_BoolShares_cross(
    const size_t N, const size_t a, const size_t b,
    const bool* x, const bool* y, bool* const z) {
  int sent_bytes = 0;

  check_BoolTriples(2 * N * a * b);

  bool* x_ext = new bool[N * a * b];
  bool* y_ext = new bool[N * a * b];
  bool* de = new bool[2 * N * a * b];
  bool* de_other = new bool[2 * N * a * b];

  multiply_BoolShares_cross_setup(N, a, b, x, y, z, x_ext, y_ext, de);

  sent_bytes += send_bool_batch(serverfd, de, 2 * N * a * b);
  recv_bool_batch(serverfd, de_other, 2 * N * a * b);

  multiply_BoolShares_finish(N * a * b, x_ext, y_ext, z, de, de_other);

  delete[] x_ext;
  delete[] y_ext;

  return sent_bytes;
}

// c_{i+1} = c_i xor ((x_i xor c_i) and (y_i xor c_i))
// output z_i = x_i xor y_i xor c_i
// 1 Bool Trip per bit
int CorrelatedStore::add_BinaryShares(
    const size_t N, const size_t* const num_bits,
    const bool* const * const x, const bool* const * const y,
    bool* const * const z, bool* const carry) {
  int sent_bytes = 0;

  memset(carry, false, N);
  bool* const xi = new bool[N];
  bool* const yi = new bool[N];

  size_t max_bits = 0;
  size_t total_bits = 0;
  for (unsigned int i = 0; i < N; i++) {
    max_bits = (num_bits[i] > max_bits ? num_bits[i] : max_bits);
    total_bits += num_bits[i];
  }

  check_BoolTriples(total_bits);

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

    bool* const new_carry = new bool[idx];
    sent_bytes += multiply_BoolShares(idx, xi, yi, new_carry);

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

  return sent_bytes;
}

// Pure OT, no correlated used
int CorrelatedStore::multiply_BoolArith(
    const size_t N, const size_t B, const bool* const b, const fmpz_t* const x,
    fmpz_t* const z, fmpz_t* const z_inv, const bool* const valid
) {
  int sent_bytes = 0;
  /* b0, b1 bits. x0, x1 values
  want z0 + z1 = (b0 ^ b1) * (x0 + x1)
  Send (r, r + x), flip order on b0, select with b1.
  So receive r if b0 = b1, r + x if b0 != b1. => r + (b0 ^ b1) x
  */
  uint64_t* const data0 = new uint64_t[N * B];
  uint64_t* const data1 = new uint64_t[N * B];
  uint64_t* const data0_inv = z_inv ? new uint64_t[N * B] : nullptr;
  uint64_t* const data1_inv = z_inv ? new uint64_t[N * B] : nullptr;
  fmpz_t r; fmpz_init(r);
  fmpz_t rx; fmpz_init(rx);
  for (unsigned int i = 0; i < N; i++) {
    if (valid && !valid[i]) continue;
    for (unsigned int j = 0; j < B; j++) {
      const unsigned int idx = i * B + j;

      fmpz_randm(r, seed, Int_Modulus);
      // fmpz_set_ui(r, 100);
      fmpz_mod_sub(z[idx], z[idx], r, mod_ctx);
      fmpz_mod_add(rx, r, x[idx], mod_ctx);
      // (r, r+x), swap if b
      data0[idx] = fmpz_get_ui(b[idx] ? rx : r);
      data1[idx] = fmpz_get_ui(b[idx] ? r : rx);

      if (z_inv) {
        fmpz_randm(r, seed, Int_Modulus);
        // fmpz_set_ui(r, 100);
        fmpz_mod_sub(z_inv[idx], z_inv[idx], r, mod_ctx);
        fmpz_mod_add(rx, r, x[idx], mod_ctx);
        // (r+x, r), swap if b
        data0_inv[idx] = fmpz_get_ui(b[idx] ? r : rx);
        data1_inv[idx] = fmpz_get_ui(b[idx] ? rx : r);
      }
    }
  }
  fmpz_clear(r);
  fmpz_clear(rx);

  uint64_t* const received = new uint64_t[N*B];
  uint64_t* const received_inv = z_inv ? new uint64_t[N*B] : nullptr;
  // Fork stuff ignored for now.
  // std::cout << "\tmultBoolARith OT size: " << (N*B) << std::endl;
  if (server_num == 0) {
    sent_bytes += ot0->send(data0, data1, N*B, data0_inv, data1_inv);
    sent_bytes += ot1->recv(received, b, N*B, received_inv);
  } else {
    sent_bytes += ot0->recv(received, b, N*B, received_inv);
    sent_bytes += ot1->send(data0, data1, N*B, data0_inv, data1_inv);
  }
  delete[] data0;
  delete[] data1;
  if (z_inv) {
    delete[] data0_inv;
    delete[] data1_inv;
  }

  for (unsigned int i = 0; i < N*B; i++) {
    fmpz_mod_add_ui(z[i], z[i], received[i], mod_ctx);
  }
  delete[] received;
  if (z_inv) {
    for (unsigned int i = 0; i < N*B; i++) {
      fmpz_mod_add_ui(z_inv[i], z_inv[i], received_inv[i], mod_ctx);
    }
    delete[] received_inv;
  }

  return sent_bytes;
}

int CorrelatedStore::multiply_BoolArithFlat(
    const size_t N, const size_t B, const bool* const b_flat, const fmpz_t* const x,
    fmpz_t* const z, fmpz_t* const z_inv, const bool* const valid
) {
  bool* b = new bool[N * B];
  for (unsigned int j = 0; j < B; j++)
    for (unsigned int i = 0; i < N; i++)
      b[i * B + j] = b_flat[j];
  int sent_bytes = multiply_BoolArith(N, B, b, x, z, z_inv, valid);
  delete[] b;
  return sent_bytes;
}

int CorrelatedStore::heavy_convert(
    const size_t N, const size_t b,
    const bool* const x, const bool* const y,
    const bool* const valid,
    fmpz_t* const bucket0, fmpz_t* const bucket1) {
  int sent_bytes = 0;

  check_DaBits(N * b);

  // Step 1: convert y to arith shares
  fmpz_t* y_p; new_fmpz_array(&y_p, N * b);
  sent_bytes += b2a_single(N * b, y, y_p);

  // Step 2: OT setup
  // z = 1 - 2y, as [z] = servernum - 2[y]
  // then (x ^ x')(z + z')
  fmpz_t* z; new_fmpz_array(&z, N * b);

  // [ (r0, r1+z), (r0 + z, r1) ] based on x
  // [ (r0 + xz, r1 + !x z), (r0 + !x z, r1 + x z)]
  // [ (0, 0_1), (1, 1_1)]
  for (unsigned int i = 0; i < N; i++) {
    if (!valid[i]) continue;

    for (unsigned int j = 0; j < b; j++) {
      const size_t idx = i * b + j;

      // Build z = 1 - 2y, as [z] = servernum - 2[y]
      fmpz_set_si(z[idx], server_num);
      fmpz_mod_submul_ui(z[idx], y_p[idx], 2, mod_ctx);
    }
  }
  clear_fmpz_array(y_p, N * b);

  fmpz_t* buff0; new_fmpz_array(&buff0, N * b);
  fmpz_t* buff1; new_fmpz_array(&buff1, N * b);
  sent_bytes += multiply_BoolArith(N, b, x, z, buff1, buff0, valid);
  clear_fmpz_array(z, N * b);
  for (unsigned int i = 0; i < N; i++) {
    for (unsigned int j = 0; j < b; j++) {
      const size_t idx = i * b + j;
      fmpz_mod_add(bucket0[j], bucket0[j], buff0[idx], mod_ctx);
      fmpz_mod_add(bucket1[j], bucket1[j], buff1[idx], mod_ctx);
    }
  }
  clear_fmpz_array(buff0, N * b);
  clear_fmpz_array(buff1, N * b);
  return sent_bytes;
}

const unsigned int HEAVY_BATCH_SIZE_BASE = 1200000;

int CorrelatedStore::heavy_convert_mask(
      const size_t N, const size_t Q, const size_t M, const size_t D,
      const bool* const x, const fmpz_t* const y_p, const bool* const mask,
      const bool* const valid, fmpz_t* const bucket0, fmpz_t* const bucket1) {
  int sent_bytes = 0;

  const size_t heavy_batch_size = HEAVY_BATCH_SIZE_BASE / (Q * M * D);

  if (N > heavy_batch_size) {
    std::cout << "heavy convert total: " << N << ", heavy batch size: " << heavy_batch_size << std::endl;
    size_t num_processed = 0;
    size_t batch_idx = 0;
    size_t num_batches = ceil(1.0 * N / heavy_batch_size);
    while (num_processed < N) {
      const size_t num_remaining = N - num_processed;
      const size_t this_batch_size = num_remaining < heavy_batch_size ? num_remaining : heavy_batch_size;
      batch_idx++;
      std::cout << "Starting heavy batch: " << batch_idx << " / " << num_batches << std::endl;

      sent_bytes += heavy_convert_mask(this_batch_size, Q, M, D,
        &x[num_processed * Q * D],
        &y_p[num_processed * Q * D],
        &mask[num_processed * Q * M],
        &valid[num_processed],
        bucket0, bucket1
      );

      num_processed += this_batch_size;
    }

    return sent_bytes;
  }

  auto start = clock_start();

  fmpz_t* z_base; new_fmpz_array(&z_base, N * Q * M * D);
  bool* x_extended = new bool[N * Q * M * D];
  bool* mask_extended = new bool[N * Q * M * D];

  fmpz_t z; fmpz_init(z);

  for (unsigned int n = 0; n < N; n++) {
    if (!valid[n]) continue;
    for (unsigned int q = 0; q < Q; q++) {
      const size_t q_idx = n * Q + q;
      for (unsigned int d = 0; d < D; d++) {
        const size_t xy_idx = q_idx * D + d;
        // Build z = 1 - 2y, as [z] = servernum - 2[y]
        fmpz_set_si(z, server_num);
        fmpz_mod_submul_ui(z, y_p[xy_idx], 2, mod_ctx);

        for (unsigned int m = 0; m < M; m++) {
          const size_t idx = (q_idx * M + m) * D + d;
          fmpz_set(z_base[idx], z);
          x_extended[idx] = x[xy_idx];
        }
      }
      for (unsigned int m = 0; m < M; m++) {
        const size_t mask_idx = q_idx * M + m;
        memset(&mask_extended[mask_idx * D], mask[mask_idx], D);
      }
    }
  }
  fmpz_clear(z);

  // Round 1: mask * z
  // Note: Can't straight multiply masks, since this does tricks for x * z
  fmpz_t* z_masked; new_fmpz_array(&z_masked, N * Q * M * D);
  sent_bytes += multiply_BoolArith(N, Q * M * D, mask_extended, z_base, z_masked, nullptr, valid);
  delete[] mask_extended;
  clear_fmpz_array(z_base, N * Q * M * D);
  std::cout << "  mask mul time: " << sec_from(start) << "\n"; start = clock_start();

  // Round 2: x * z
  // Uses tricks for second buff, to implicitly map x to (x, 1-x)
  // By assuming second one is always "inverted".
  // Would require extra to allow for "both 0" case with mask
  fmpz_t* buff0; new_fmpz_array(&buff0, N * Q * M * D);
  fmpz_t* buff1; new_fmpz_array(&buff1, N * Q * M * D);
  sent_bytes += multiply_BoolArith(N, Q * M * D, x_extended, z_masked, buff1, buff0, valid);
  delete[] x_extended;
  clear_fmpz_array(z_masked, N * Q * M * D);
  std::cout << "  xz mul time: " << sec_from(start) << "\n"; start = clock_start();

  // for (unsigned int i = 0; i < 1 + 0 * N * Q * M * D; i++) {
  //   std::cout << "x_ext[" << i << "]_" << server_num << " = " << x_extended[i] << std::endl;
  //   std::cout << "mask_ext[" << i << "]_" << server_num << " = " << mask_extended[i] << std::endl;
  //   std::cout << "z_base[" << i << "]_" << server_num << " = " << get_fsigned(z_base[i], Int_Modulus) << std::endl;
  //   std::cout << "z_masked[" << i << "]_" << server_num << " = " << get_fsigned(z_masked[i], Int_Modulus) << std::endl;
  //   std::cout << "buff0[" << i << "]_" << server_num << " = " << get_fsigned(buff0[i], Int_Modulus) << std::endl;
  //   std::cout << "buff1[" << i << "]_" << server_num << " = " << get_fsigned(buff1[i], Int_Modulus) << std::endl;
  // }

  for (unsigned int bucket_idx = 0; bucket_idx < Q * M * D; bucket_idx++) {
    for (unsigned int n = 0; n < N; n++) {
      const size_t idx = n * (Q * M * D) + bucket_idx;
      fmpz_mod_add(bucket0[bucket_idx], bucket0[bucket_idx], buff0[idx], mod_ctx);
      fmpz_mod_add(bucket1[bucket_idx], bucket1[bucket_idx], buff1[idx], mod_ctx);
    }
  }
  clear_fmpz_array(buff0, N * Q * M * D);
  clear_fmpz_array(buff1, N * Q * M * D);

  return sent_bytes;
}


// WARNING: Just does x in clear. Early version for debug/test
// [x] > c for known c, shares [x]
// Treats > N/2 as negative
bool* CorrelatedStore::cmp_c_clear(const size_t N,
                                   const fmpz_t* const x,
                                   const fmpz_t* const c) {
  bool* const ans = new bool[N];

  fmpz_t half; fmpz_init(half); fmpz_cdiv_q_ui(half, Int_Modulus, 2);

  // Done in clear
  if (server_num == 0) {
    fmpz_t* x2; new_fmpz_array(&x2, N);
    fmpz_t* c2; new_fmpz_array(&c2, N);
    recv_fmpz_batch(serverfd, x2, N);
    for (unsigned int i = 0; i < N; i++) {
      fmpz_mod_add(x2[i], x2[i], x[i], mod_ctx);
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

// Just (x < y) = (x - y < 0) = is_neg([x - y])
int CorrelatedStore::cmp(const size_t N,
                         const fmpz_t* const x, const fmpz_t* const y,
                         fmpz_t* const ans) {
  int sent_bytes = 0;
  check_DaBits(4 * N * nbits_mod);
  check_Triples(13 * N * nbits_mod);

  fmpz_t* diff; new_fmpz_array(&diff, N);

  for (unsigned int i = 0; i < N; i++) {
    fmpz_mod_sub(diff[i], x[i], y[i], mod_ctx);
  }

  sent_bytes += is_negative(N, diff, ans);
  clear_fmpz_array(diff, N);
  return sent_bytes;
}

// |x| = (1 - 2[x < 0]) * x
// if hashes known, sign might reveal some info about non-heavy bucket, so reveal in clear is not fully secure
int CorrelatedStore::abs(const size_t N, const fmpz_t* const x,
                         fmpz_t* const abs_x) {
  int sent_bytes = 0;
  check_DaBits(4 * N * nbits_mod);
  check_Triples((13 + 1) * N * nbits_mod);

  fmpz_t* is_neg; new_fmpz_array(&is_neg, N);
  sent_bytes += is_negative(N, x, is_neg);
  for (unsigned int i = 0; i < N; i++) {
    fmpz_mod_mul_si(is_neg[i], is_neg[i], -2, mod_ctx);
    if (server_num == 1)
      fmpz_mod_add_ui(is_neg[i], is_neg[i], 1, mod_ctx);
  }
  sent_bytes += multiply_ArithmeticShares(N, x, is_neg, abs_x);
  clear_fmpz_array(is_neg, N);
  return sent_bytes;
}

int CorrelatedStore::abs_cmp(const size_t N,
                             const fmpz_t* const x, const fmpz_t* const y,
                             fmpz_t* const ans) {
  int sent_bytes = 0;
  check_DaBits(3 * 4 * N * nbits_mod);  // 1 per abs, and 1 for cmp
  check_Triples((13 + 14 + 14) * N * nbits_mod);  // 13 for cmp, 14 for abs

  // Do abs x and y in one merged set
  fmpz_t* merge; new_fmpz_array(&merge, 2*N);
  for (unsigned int i = 0; i < N; i++) {
    fmpz_set(merge[i], x[i]);
    fmpz_set(merge[i+N], y[i]);
  }
  fmpz_t* merge_abs; new_fmpz_array(&merge_abs, 2*N);
  sent_bytes += abs(2*N, merge, merge_abs);
  clear_fmpz_array(merge, 2*N);

  // split out, and compare
  fmpz_t* x2; new_fmpz_array(&x2, N);
  fmpz_t* y2; new_fmpz_array(&y2, N);
  for (unsigned int i = 0; i < N; i++) {
    fmpz_set(x2[i], merge_abs[i]);
    fmpz_set(y2[i], merge_abs[i+N]);
  }
  clear_fmpz_array(merge_abs, N);
  sent_bytes += cmp(N, x2, y2, ans);
  clear_fmpz_array(x2, N);
  clear_fmpz_array(y2, N);
  return sent_bytes;
}

// x, y are Nxb shares of N total b-bit numbers.
// I.e. x[i,j] is additive share of bit j of number i.
int CorrelatedStore::cmp_bit(const size_t N, const size_t b,
                             const fmpz_t* const x, const fmpz_t* const y,
                             fmpz_t* const ans) {
  int sent_bytes = 0;
  size_t idx;  // for convenience

  // Total mults: ~3x (N*b)
  check_Triples(3 * N * b);

  // [c] = [x ^ y] = [x] + [y] - 2[xy]
  fmpz_t* c; new_fmpz_array(&c, N * b);
  sent_bytes += multiply_ArithmeticShares(N * b, x, y, c);
  for (unsigned int i = 0; i < N * b; i++) {
    fmpz_mod_mul_si(c[i], c[i], -2, mod_ctx);
    fmpz_mod_add(c[i], c[i], x[i], mod_ctx);
    fmpz_mod_add(c[i], c[i], y[i], mod_ctx);
  }

  // [di] = OR([cj]) from i+1 to b
  // = [ci] OR [d(i+1)], with d_b = cb
  // a OR b = 1-(1-a)(1-b) = a + b - ab
  // TODO: Currently doing b-round version.
  //  Can be constant, but lazy fine for now.
  //  wait, can we just fold it into one large multiply?
  //  Set ci, di1. mul = ci * di1. Then set d
  fmpz_t* d; new_fmpz_array(&d, N * b);
  // Start: db = cb
  for (unsigned int i = 0; i < N; i++) {
    idx = i * b + (b-1);  // [i, b-1]
    fmpz_set(d[idx], c[idx]);
  }
  fmpz_t* ci; new_fmpz_array(&ci, N);     // [ci]
  fmpz_t* di1; new_fmpz_array(&di1, N);   // [d(i+1)]
  fmpz_t* mul; new_fmpz_array(&mul, N);
  for (int j = b-2; j >= 0; j--) {
    for (unsigned int i = 0; i < N; i++) {
      idx = i * b + j;
      fmpz_set(ci[i], c[idx]);
      fmpz_set(di1[i], d[idx + 1]);
    }
    sent_bytes += multiply_ArithmeticShares(N, ci, di1, mul);
    for (unsigned int i = 0; i < N; i++) {
      idx = i * b + j;
      fmpz_mod_add(d[idx], ci[i], di1[i], mod_ctx);
      fmpz_mod_sub(d[idx], d[idx], mul[i], mod_ctx);
    }
  }
  clear_fmpz_array(c, N * b);
  clear_fmpz_array(ci, N);
  clear_fmpz_array(di1, N);
  clear_fmpz_array(mul, N);

  // ei = di - d(i+1), with eb = db
  fmpz_t* e; new_fmpz_array(&e, N * b);
  for (unsigned int i = 0; i < N; i++) {
    idx = i * b + (b-1);
    fmpz_set(e[idx], d[idx]);
    for (unsigned int j = 0; j < b - 1; j++) {
      idx = i * b + j;
      fmpz_mod_sub(e[idx], d[idx], d[idx + 1], mod_ctx);
    }
  }
  clear_fmpz_array(d, N * b);

  // [x < y] = sum ei * yi
  fmpz_t* ey; new_fmpz_array(&ey, N * b);
  sent_bytes += multiply_ArithmeticShares(N * b, e, y, ey);
  clear_fmpz_array(e, N * b);
  for (unsigned int i = 0; i < N; i++) {
    fmpz_zero(ans[i]);
    for (unsigned int j = 0; j < b; j++) {
      fmpz_mod_add(ans[i], ans[i], ey[i * b + j], mod_ctx);
    }
  }
  clear_fmpz_array(ey, N * b);

  return sent_bytes;
}

// Make each bit in parallel: [ri in {0,1}]p
// check if [r < p] bitwise, retry if not
// success odds are p/2^b, worst case ~1/2 failure.
int CorrelatedStore::gen_rand_bitshare(const size_t N,
                                       fmpz_t* const r, fmpz_t* const rB) {
  int sent_bytes = 0;
  const size_t b = nbits_mod;

  // Assumed tries to succeed. (worst case 1/2 fail, so 2 avg, overestimate)
  const size_t avg_tries = 4;
  check_DaBits(avg_tries * N * b);       // for gen
  check_Triples(avg_tries * 3 * N * b);  // for cmp

  // p bitwise shared, for [r < p]. All bits set on server 0, since public.
  fmpz_t* pB; new_fmpz_array(&pB, N * b);  // Zero by default
  if (server_num == 0)
    for (unsigned int j = 0; j < b; j++)
      if (fmpz_tstbit(Int_Modulus, j))
        for (unsigned int i = 0; i < N; i++)
          fmpz_set_si(pB[i * b + j], 1);

  // for validation checking
  fmpz_t* rB_tocheck; new_fmpz_array(&rB_tocheck, N * b);
  size_t* rB_idx = new size_t[N];
  fmpz_t* r_lt_p_other; new_fmpz_array(&r_lt_p_other, N);

  // Track which are already set
  bool* const valid = new bool[N];
  memset(valid, false, N * sizeof(bool));
  size_t num_invalid;
  [[maybe_unused]] size_t log_num_invalid = 0;  // Just for logging

  while (true) {
    num_invalid = 0;

    check_DaBits(N * b);

    // Compute new (if not valid)
    for (unsigned int i = 0; i < N; i++) {
      if (valid[i])
        continue;

      fmpz_zero(r[i]);
      for (unsigned int j = 0; j < b; j++) {
        // Just need 2 numbers summing to 0 or 1, so bp. b2 not needed.
        const DaBit* const dabit = get_DaBit();
        fmpz_set(rB[i * b + j], dabit->bp);
        fmpz_mod_addmul_ui(r[i], dabit->bp, 1ULL << j, mod_ctx);
        // consume the dabit
        delete dabit;

        // add to rB_tocheck
        // std::cout << "rB_tocheck[" << num_invalid * b + j << "] := rB[" << i*b+j << "]" << std::endl;
        fmpz_set(rB_tocheck[num_invalid * b + j], rB[i * b + j]);
      }
      rB_idx[num_invalid] = i;

      num_invalid++;
    }

    log_num_invalid += num_invalid;

    // std::cout << "num invalid: " << num_invalid << " / " << N << std::endl;
    if (num_invalid == 0) {
      break;
    }

    // Check [r < p], retry if not, sets "valid" where [r < p]

    // It's fine if arrays are larger, extras get ignored.
    fmpz_t* r_lt_p; new_fmpz_array(&r_lt_p, num_invalid);
    sent_bytes += cmp_bit(num_invalid, b, rB_tocheck, pB, r_lt_p);
    // Get r_lt_p in clear
    sent_bytes += reveal_fmpz_batch(serverfd, r_lt_p, num_invalid);
    for (unsigned int i = 0; i < num_invalid; i++) {
      valid[rB_idx[i]] = fmpz_is_one(r_lt_p[i]);
      // std::cout << "check " << i << ": valid[" << rB_idx[i] << "] = " << valid[rB_idx[i]] << std::endl;
    }
    clear_fmpz_array(r_lt_p, num_invalid);
  }

  clear_fmpz_array(rB_tocheck, N * b);
  clear_fmpz_array(pB, N * b);
  clear_fmpz_array(r_lt_p_other, N);
  delete[] rB_idx;
  delete[] valid;

  // std::cout << "total sub-iterations: " << log_num_invalid << ", vs expected: " << avg_tries * N << std::endl;

  return sent_bytes;
}

// Since we are masking with random r, this should have max b (2^b > p)
int CorrelatedStore::LSB(const size_t N, const fmpz_t* const x,
                         fmpz_t* const x0) {
  int sent_bytes = 0;
  const size_t b = nbits_mod;

  check_DaBits(4 * N * b);   // for gen_rand
  check_Triples((4*3 + 1) * N * b);  // 12 for gen_rand

  // 1: Random bitwise shared r, true [r]p, bitwise [rB]
  fmpz_t* r; new_fmpz_array(&r, N);
  fmpz_t* rB; new_fmpz_array(&rB, N * b);
  sent_bytes += gen_rand_bitshare(N, r, rB);

  // 2: Compute [c]p = [x]p + [r]p
  fmpz_t* c; new_fmpz_array(&c, N);
  for (unsigned int i = 0; i < N; i++) {
    fmpz_mod_add(c[i], x[i], r[i], mod_ctx);
  }
  clear_fmpz_array(r, N);
  // 2.1: Reveal c = x + r
  sent_bytes += reveal_fmpz_batch(serverfd, c, N);

  // 3: if wraparound (c < r), x0 = c0 ^ r0
  //    if no wraparound, then x0 = 1 - c0 ^ r0
  // [x0] = [c < r] + (c0 ^ r0) + 2 [c < r] (c0 ^ r0)
  // 3.1: c0 ^ r0 = ([r0] if c0 = 0, 1 - [r0] if c0 = 1
  fmpz_t* cr; new_fmpz_array(&cr, N);
  fmpz_t adj; fmpz_init_set_ui(adj, server_num);  // "1" is just once
  for (unsigned int i = 0; i < N; i++) {
    fmpz_set(cr[i], rB[i * b]);  // Just get r0 from bitwise shares
    if (fmpz_tstbit(c[i], 0) == 1) {  // 1 - cr
      fmpz_mod_sub(cr[i], adj, cr[i], mod_ctx);
    }
  }
  fmpz_clear(adj);
  // 3.2: Get [c < r], as bitwise
  // Since c is in clear, server 0 sets [c] = c, 1 sets [c] = 0
  fmpz_t* cB; new_fmpz_array(&cB, N * b);
  if (server_num == 0) {
    for (unsigned int i = 0; i < N; i++) {
      for (unsigned int j = 0; j < b; j++) {
        fmpz_set_ui(cB[i * b + j], fmpz_tstbit(c[i], j));
      }
    }
  } else {
    for (unsigned int i = 0; i < N * b; i++) {
      fmpz_zero(cB[i]);
    }
  }
  fmpz_t* cmp; new_fmpz_array(&cmp, N);
  sent_bytes += cmp_bit(N, b, cB, rB, cmp);
  clear_fmpz_array(cB, N*b);
  clear_fmpz_array(rB, N*b);
  clear_fmpz_array(c, N);
  // 4: [x0] = [c < r] + (c0 ^ r0) - 2 [c < r] (c0 ^ r0)
  // 4.1: mul = [c<r] * [c0 ^ r0]
  fmpz_t* mul; new_fmpz_array(&mul, N);
  sent_bytes += multiply_ArithmeticShares(N, cmp, cr, mul);
  // 4.2: Final eval: cmp + cr - 2 cmp*cr
  for (unsigned int i = 0; i < N; i++) {
    fmpz_mod_add(x0[i], cmp[i], cr[i], mod_ctx);
    fmpz_mod_submul_ui(x0[i], mul[i], 2, mod_ctx);
  }

  clear_fmpz_array(cmp, N);
  clear_fmpz_array(cr, N);
  clear_fmpz_array(mul, N);

  return sent_bytes;
}

int CorrelatedStore::is_negative(const size_t N,
                                 const fmpz_t* const x, fmpz_t* ans) {
  int sent_bytes = 0;
  check_DaBits(4 * N * nbits_mod);
  check_Triples(13 * N * nbits_mod);

  // [(2a)_0]
  fmpz_t* inp; new_fmpz_array(&inp, N);
  for (unsigned int i = 0; i < N; i++) {
    fmpz_mod_mul_ui(inp[i], x[i], 2, mod_ctx);
  }
  sent_bytes += LSB(N, inp, ans);
  clear_fmpz_array(inp, N);
  // Code for (1 - [(2a)_0]), aka 1 if positive, 0 if negative
  // for (unsigned int i = 0; i < N; i++) {
  //   fmpz_mod_neg(ans[i], ans[i], mod_ctx);
  //   if (server_num == 0)
  //     fmpz_mod_add_ui(ans[i], ans[i], 1, mod_ctx);
  // }
  return sent_bytes;
}



////////////////

int PrecomputeStore::check_DaBits(const size_t n) {
  if (dabit_store.size() < n)
    return add_DaBits(n - dabit_store.size());
  return 0;
}

void PrecomputeStore::check_BoolTriples(const size_t n) {
  if (btriple_store.size() < n) add_BoolTriples(n - btriple_store.size());
}

void PrecomputeStore::check_Triples(const size_t n) {
  if (atriple_store.size() < n) add_Triples(n - atriple_store.size());
}

int PrecomputeStore::add_DaBits(const size_t n) {
  auto start = clock_start();
  int sent_bytes = 0;
  const size_t num_to_make = (n > batch_size ? n : batch_size);
  // const size_t num_to_make = n;  // Use this to make "end to end" easier to benchmark
  std::cout << "adding dabits: " << num_to_make << std::endl;
  DaBit** const dabits = new DaBit*[num_to_make];

  if (lazy) {
    sent_bytes += gen_DaBits_lazy(num_to_make, dabits);
  } else {
    sent_bytes += gen_DaBits(num_to_make, dabits);
  }
  for (unsigned int i = 0; i < num_to_make; i++)
    dabit_store.push(dabits[i]);
  delete[] dabits;
  std::cout << "add_DaBits timing : " << sec_from(start) << std::endl;
  return sent_bytes;
}

void PrecomputeStore::add_Triples(const size_t n) {
  auto start = clock_start();
  const size_t num_to_make = (n > batch_size ? n : batch_size);
  std::cout << "adding triples: " << num_to_make << std::endl;
  // if (triple_gen) {  // not null pointer
  //   std::vector<BeaverTriple*> new_triples = triple_gen->generateTriples(num_to_make);
  //   for (unsigned int i = 0; i < num_to_make; i++)
  //     atriple_store.push(new_triples[i]);
  // } else {
  std::cout << "Using lazy beaver triples" << std::endl;
  for (unsigned int i = 0; i < num_to_make; i++) {
    const BeaverTriple* const triple = generate_beaver_triple_lazy(serverfd, server_num);
    atriple_store.push(triple);
    // }
  }
  std::cout << "add_Triples timing : " << sec_from(start) << std::endl;
}

void PrecomputeStore::add_BoolTriples(const size_t n) {
  auto start = clock_start();
  const size_t num_to_make = (n > batch_size ? n : batch_size);
  std::cout << "adding booltriples: " << num_to_make << std::endl;
  std::queue<const BooleanBeaverTriple*> new_triples = gen_boolean_beaver_triples(server_num, num_to_make, ot0, ot1);
  for (unsigned int i = 0; i < num_to_make; i++) {
    btriple_store.push(new_triples.front());
    new_triples.pop();
  }
  std::cout << "add_BoolTriples timing : " << sec_from(start) << std::endl;
}

void PrecomputeStore::print_Sizes() const {
  std::cout << "Current store sizes:\n";
  std::cout << " Dabits: " << dabit_store.size() << std::endl;
  std::cout << " Bool  Triples: " << btriple_store.size() << std::endl;
  std::cout << " Arith Triples: " << atriple_store.size() << std::endl;
}

void PrecomputeStore::maybe_Update() {
  auto start = clock_start();

  // Extra make factor, to overkill for convenience
  const size_t extra_make_factor = 2;
  // Targets.
  const size_t da_target = batch_size;
  const size_t atrip_target = batch_size;
  // Estimating factor of 32. based on num bits, specific use, etc.
  const size_t btrip_target = batch_size * 32;

  if (btriple_store.size() < btrip_target)
    add_BoolTriples(btrip_target * extra_make_factor);

  if (dabit_store.size() < da_target)
    add_DaBits(da_target * extra_make_factor);

  if (atriple_store.size() < atrip_target)
    add_Triples(atrip_target * extra_make_factor);

  // print_Sizes();

  std::cout << "precompute timing : " << sec_from(start) << std::endl;
}

// Use b2A via OT on random bit
// Nearly COT, except delta is changing
// random choice and random base, but also random delta matters
// TODO: Work on larger values. Currently assumes mod is uint64_t.
int PrecomputeStore::gen_DaBits(const size_t N, DaBit** const dabit) {
  int sent_bytes = 0;
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
    sent_bytes += ot0->send(b0, b1, N);
    delete[] b0;
    delete[] b1;
  } else {
    sent_bytes += ot0->recv(x, b, N);
  }
  for (unsigned int i = 0; i < N; i++) {
    uint64_t bp = (b[i] + 2 * (mod - x[i])) % mod;
    fmpz_set_ui(dabit[i]->bp, bp);
  }

  delete[] b;
  delete[] x;

  return sent_bytes;
}

int PrecomputeStore::gen_DaBits_lazy(const size_t N, DaBit** const dabit) {
  int sent_bytes = 0;
  for (unsigned int i = 0; i < N; i++)
    dabit[i] = new DaBit();
  if (server_num == 0) {
    DaBit** const dabit_other = new DaBit*[N];
    for (unsigned int i = 0; i < N; i++) {
      dabit_other[i] = new DaBit();
      makeLocalDaBit(dabit[i], dabit_other[i]);
    }
    sent_bytes += send_DaBit_batch(serverfd, dabit_other, N);
    for (unsigned int i = 0; i < N; i++)
      delete dabit_other[i];
    delete[] dabit_other;
  } else {
    recv_DaBit_batch(serverfd, dabit, N);
  }
  return sent_bytes;
}



////////////////


int OTCorrelatedStore::b2a_multi(
    const size_t N, const size_t* const num_bits,
    const fmpz_t* const x, fmpz_t* const xp) {
  return b2a_ot(1, N, num_bits, x, xp, fmpz_get_ui(Int_Modulus));
}

int OTCorrelatedStore::b2a_single(
    const size_t N, const bool* const x, fmpz_t* const xp) {
  size_t num_bits[N];
  fmpz_t* x2; new_fmpz_array(&x2, N);
  for (unsigned int i = 0; i < N; i++) {
    fmpz_set_ui(x2[i], (int) x[i]);
    num_bits[i] = 1;
  }
  int ret = b2a_multi(N, num_bits, x2, xp);
  clear_fmpz_array(x2, N);
  return ret;
}

// Using intsum_ot, multiple bits
int OTCorrelatedStore::b2a_ot(
    const size_t num_shares, const size_t num_values, const size_t* const num_bits,
    const fmpz_t* const x, fmpz_t* const xp_out, const size_t mod) {
  int sent_bytes = 0;

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

  // TODO: incorperate sent_bytes into OT
  if (server_num == 0) {
    xp = intsum_ot_sender(ot0, x2, valid, num_bits, num_shares, num_values, mod);
  } else {
    xp = intsum_ot_receiver(ot0, x2, num_bits, num_shares, num_values, mod);
  }

  // for consistency, flatten and fmpz_t
  for (unsigned int i = 0; i < num_shares; i++) {
    for (unsigned int j = 0; j < num_values; j++) {
      fmpz_set_ui(xp_out[i * num_values + j], xp[i][j]);
    }
    delete[] x2[i];
    delete[] xp[i];
  }
  delete[] x2;
  delete[] valid;
  delete[] xp;

  return sent_bytes;
}
