#undef NDEBUG
#include <iostream>

#include "utils_test_connect.h"
#include "../constants.h"
#include "../correlated.h"
#include "../fmpz_utils.h"
#include "../net_share.h"


const size_t batch_size = 100; // flexible
const int lazy = 0;

const bool use_ot_version = false;
const bool use_large = false;

/*
4 bits, for each bucket/hash combo
8 inputs per, for possible randomness (0-bucket, 2 xor bits)
Tries all combinations of shares

bit = 0, all +1 into bucket 0
bit = 1, all +1 into bucket 1
bit = 2, all -1 into bucket 0
bit = 3, all -1 into bucket 1
*/

void test_HeavyConvert(const int server_num, const int serverfd,
                       CorrelatedStore* store) {
  const size_t nbits = 16;
  const size_t N = use_large ? 100000 : 100;
  const size_t n = nbits * N;

  bool* const x = new bool[n];
  bool* const y = new bool[n];
  fmpz_t* bucket0; new_fmpz_array(&bucket0, nbits);
  fmpz_t* bucket1; new_fmpz_array(&bucket1, nbits);
  bool* const valid = new bool[N]; memset(valid, true, N);

  int expected0[nbits]; memset(expected0, 0, nbits * sizeof(int));
  int expected1[nbits]; memset(expected1, 0, nbits * sizeof(int));

  // Setup
  for (unsigned int j = 0; j < nbits; j++) {
    const bool bucket = j % 2;
    const bool hash  = (j >> 1) % 2;
    // if (server_num == 0) {
    //   std::cout << "bit " << j << " = bucket " << bucket << ", hash " << hash << std::endl;
    // }

    (bucket ? expected1 : expected0)[j] = (hash ? -(int)N : N);

    for (unsigned int i = 0; i < N; i++) {
      const size_t idx = i * nbits + j;

      // share 0
      x[idx] = i % 2;
      y[idx] = (i >> 1) % 2;
      // share 1
      if (server_num == 1) {
        // std::cout << "x[" << idx << "] = " << x[idx] << ", " << (x[idx] ^ bucket) << std::endl;
        // std::cout << "y[" << idx << "] = " << y[idx] << ", " << (y[idx] ^ hash) << std::endl;
        x[idx] ^= bucket;
        y[idx] ^= hash;
      }
    }
  }

  if (use_ot_version) {
    store->check_DaBits(N * nbits);
  } else {
    store->check_BoolTriples(2 * N * nbits);
    store->check_DaBits(3 * N * nbits);
  }

  auto start = clock_start();
  std::cout << "clock start" << std::endl;
  int sent_bytes = 0;
  if (use_ot_version) {
    sent_bytes = store->heavy_convert_ot(N, nbits, x, y, valid, bucket0, bucket1);
  } else {
    sent_bytes = store->heavy_convert(N, nbits, x, y, valid, bucket0, bucket1);
  }
  std::cout << "heavy convert timing : " << sec_from(start) << std::endl;
  std::cout << "sent_bytes = " << sent_bytes << std::endl;

  // recombine / test
  if (server_num == 0) {
    fmpz_t* bucket0_other; new_fmpz_array(&bucket0_other, nbits);
    fmpz_t* bucket1_other; new_fmpz_array(&bucket1_other, nbits);
    recv_fmpz_batch(serverfd, bucket0_other, nbits);
    recv_fmpz_batch(serverfd, bucket1_other, nbits);
    fmpz_t tmp; fmpz_init(tmp);
    for (unsigned int j = 0; j < nbits; j++) {
      fmpz_mod_add(tmp, bucket0[j], bucket0_other[j], mod_ctx);
      // std::cout << "bucket0[" << j << "] total = " << fmpz_mod_get_signed(tmp, Int_Modulus);
      // std::cout << ", \texpected = " << expected0[j] << std::endl;
      assert(fmpz_mod_get_signed(tmp, Int_Modulus) == expected0[j]);
      fmpz_mod_add(tmp, bucket1[j], bucket1_other[j], mod_ctx);
      // std::cout << "bucket1[" << j << "] total = " << fmpz_mod_get_signed(tmp, Int_Modulus);
      // std::cout << ", \texpected = " << expected1[j] << std::endl;
      assert(fmpz_mod_get_signed(tmp, Int_Modulus) == expected1[j]);
    }
    fmpz_clear(tmp);
    clear_fmpz_array(bucket0_other, nbits);
    clear_fmpz_array(bucket1_other, nbits);
  } else {
    send_fmpz_batch(serverfd, bucket0, nbits);
    send_fmpz_batch(serverfd, bucket1, nbits);
  }

  clear_fmpz_array(bucket0, nbits);
  clear_fmpz_array(bucket1, nbits);
  delete[] x;
  delete[] y;
  delete[] valid;
}

void test_HeavyConvertMask_one(const int server_num, const int serverfd,
    CorrelatedStore* store) {
  const size_t N = use_large ? 10000 : 80;
  const size_t Q = use_large ? 13 : 3;
  const size_t M = use_large ? 11 : 5;
  const size_t D = use_large ? 8 : 4;
  // Q "parallel" runs of heavy_convert
  // N inputs, Q copies, D depth, M substreams
  // |x| = |y| = N * Q * D
  // |mask| = N * Q * M: Substream select
  // Valid size N
  // buckets size : Q * M * D
  // Order: N, Q, M, D
  //   x,y = [over N (over Q (over D))]
  //   valid = over N
  //   mask  = [over N (over Q (over M))]
  //   buckets = [over Q (over M (over D))]
  bool* const x = new bool[N * Q * D];
  bool* const y = new bool[N * Q * D];
  fmpz_t* bucket0; new_fmpz_array(&bucket0, Q * M * D);
  fmpz_t* bucket1; new_fmpz_array(&bucket1, Q * M * D);
  bool* const valid = new bool[N]; memset(valid, true, N);

  int expected0[Q * M * D]; memset(expected0, 0, Q * M * D * sizeof(int));
  int expected1[Q * M * D]; memset(expected1, 0, Q * M * D * sizeof(int));

  bool* mask = new bool[N * Q * M];

  auto mask_eval = [](int n, int m){return (n + m) % 3 == 0;};

  // Setup
  /*
  "base random"
    x, y rotate 00, 10, 01, 11
    mask at base alternates (0, 1, 0, 1)
  actual value
    x alternates along d, offset by q
    y alternates every other by d, offset by q
    mask is every 3rd
  */
  for (unsigned int n = 0; n < N; n++) {
    memset(&x[n * (Q * D)], n % 2, Q * D);
    memset(&y[n * (Q * D)], (n>>1) % 2, Q * D);
    memset(&mask[n * (Q * M)], n % 2, Q * M);
    if (server_num == 1) {
      for (unsigned int q = 0; q < Q; q++) {
        const size_t q_idx = n * Q + q;

        for (unsigned int d = 0; d < D; d++) {
          const size_t xy_idx = q_idx * D + d;
          x[xy_idx] ^= (d + q) % 2;
          y[xy_idx] ^= ((d>>1) + q) % 2;
          // std::cout << "x[" << xy_idx << "] = " << (d + q) % 2 << std::endl;
          // std::cout << "y[" << xy_idx << "] = " << ((d>>1) + q) % 2 << std::endl;
        }
        for (unsigned int m = 0; m < M; m++) {
          mask[q_idx * M + m] ^= mask_eval(n, m);
          // std::cout << "mask[" << q_idx * M + m << "] = " << mask_eval(n, m) << std::endl;
        }
      }
    }
  }

  // Setup expected.
  for (unsigned int q = 0; q < Q; q++) {
    for (unsigned int d = 0; d < D; d++) {
      // x choice
      auto e = (d + q) % 2 ? expected1 : expected0;
      for (unsigned int m = 0; m < M; m++) {
        const size_t idx = (q * M + m) * D + d;
        // 1/3 plus remainders
        e[idx] = N/3 + ((3*M - m)%3<(N%3));
        // y sign
        if (((d>>1) + q) % 2 == 1) e[idx] *= -1;
      }
    }
  }

  if (use_ot_version) {
    store->check_DaBits(N * Q * D);
  } else {
    store->check_BoolTriples(3 * N * Q * M * D);
    store->check_DaBits(4 * N * Q * M * D);
  }

  // Run
  auto start = clock_start();
  std::cout << "clock start" << std::endl;
  int64_t sent_bytes = 0;
  if (use_ot_version) {
    fmpz_t* y_p; new_fmpz_array(&y_p, N * Q * D);
    sent_bytes += store->b2a_single(N * Q * D, y, y_p);
    std::cout << "heavy mask: b2a timing " << sec_from(start) << std::endl;
    auto start2 = clock_start();
    sent_bytes += store->heavy_convert_mask(N, Q, M, D, x, y_p, mask, valid, bucket0, bucket1);
    std::cout << "heavy mask convert timing : " << sec_from(start2) << std::endl;
    std::cout << "heavy mask total timing : " << sec_from(start) << std::endl;
    clear_fmpz_array(y_p, N * Q * D);
  } else {
    sent_bytes += store->heavy_convert_mask_one(N, Q, M, D, x, y, mask, valid, bucket0, bucket1);
    std::cout << "heavy mask convert timing : " << sec_from(start) << std::endl;
  }
  std::cout << "sent_bytes = " << sent_bytes << std::endl;
  delete[] x;
  delete[] y;
  delete[] mask;
  delete[] valid;

  // Recombine / test
  if (server_num == 0) {
    fmpz_t* bucket0_other; new_fmpz_array(&bucket0_other, Q * M * D);
    fmpz_t* bucket1_other; new_fmpz_array(&bucket1_other, Q * M * D);
    recv_fmpz_batch(serverfd, bucket0_other, Q * M * D);
    recv_fmpz_batch(serverfd, bucket1_other, Q * M * D);
    fmpz_t tmp0; fmpz_init(tmp0);
    fmpz_t tmp1; fmpz_init(tmp1);
    for (unsigned int i = 0; i < Q * M * D; i++) {
      [[maybe_unused]] const size_t q = i / (M * D);
      [[maybe_unused]] const size_t m = (i / D) % M;
      [[maybe_unused]] const size_t d = i % D;

      fmpz_mod_add(tmp0, bucket0[i], bucket0_other[i], mod_ctx);
      fmpz_mod_add(tmp1, bucket1[i], bucket1_other[i], mod_ctx);

      // std::cout << "bucket[" << q << ", " << m << ", " << d << "]: ";
      // std::cout << "0 = " << fmpz_mod_get_signed(tmp0, Int_Modulus);
      // std::cout << ", expected = " << expected0[i] << ", \t";
      // std::cout << "1 = " << fmpz_mod_get_signed(tmp1, Int_Modulus);
      // std::cout << ", expected = " << expected1[i] << std::endl;

      assert(fmpz_mod_get_signed(tmp0, Int_Modulus) == expected0[i]);
      assert(fmpz_mod_get_signed(tmp1, Int_Modulus) == expected1[i]);
    }
    fmpz_clear(tmp0);
    fmpz_clear(tmp1);
    clear_fmpz_array(bucket0_other, Q * M * D);
    clear_fmpz_array(bucket1_other, Q * M * D);
  } else {
    send_fmpz_batch(serverfd, bucket0, Q * M * D);
    send_fmpz_batch(serverfd, bucket1, Q * M * D);
  }
  clear_fmpz_array(bucket0, Q * M * D);
  clear_fmpz_array(bucket1, Q * M * D);
}

void test_HeavyConvertMask_two(const int server_num, const int serverfd,
    CorrelatedStore* store) {
  const size_t N = use_large ? 10000 : 20;
  const size_t Q = use_large ? 4 : 3;
  const size_t M1 = use_large ? 9 : 2;
  const size_t M2 = use_large ? 7 : 3;
  const size_t D = use_large ? 8 : 4;

  bool* const x = new bool[N * Q * D];
  bool* const y = new bool[N * Q * D];
  bool* const valid = new bool[N]; memset(valid, true, N);

  // Note: Mask2 is currently "global" across Q.
  // Won't matter since anything it mutls against has Q. But provided this way.
  bool* mask1 = new bool[N * Q * M1];
  bool* mask2 = new bool[N * M2];

  fmpz_t* bucket0; new_fmpz_array(&bucket0, Q * M1 * M2 * D);
  fmpz_t* bucket1; new_fmpz_array(&bucket1, Q * M1 * M2 * D);
  int expected0[Q * M1 * M2 * D]; memset(expected0, 0, Q * M1 * M2 * D * sizeof(int));
  int expected1[Q * M1 * M2 * D]; memset(expected1, 0, Q * M1 * M2 * D * sizeof(int));

  // Setup
  /*
  base random:
    x, mask1 alternate player
    y, mask2 alternate every other player
  actual
    x alternates along d, offset by q
    y alternates every other d, offset by q
    m1: every third, offset by q
    m2: every third, offest by n
  */
  for (unsigned int n = 0; n < N; n++) {
    memset(&x[n * Q * D], n % 2, Q * D);
    memset(&y[n * Q * D], (n>>1) % 2, Q * D);
    memset(&mask1[n * Q * M1], n % 2, Q * M1);
    memset(&mask2[n * M2], (n>>1) % 2, M2);
    if (server_num == 1) {
      for (unsigned int q = 0; q < Q; q++) {
        const size_t q_idx = n * Q + q;

        for (unsigned int d = 0; d < D; d++) {
          x[q_idx * D + d] ^= (d + q) % 2;
          y[q_idx * D + d] ^= ((d>>1) + q) % 2;
        }

        for (unsigned int m = 0; m < M1; m++) {
          mask1[q_idx * M1 + m] ^= ( (m + q) % 3 == 0 );
        }
      }
      for (unsigned int m = 0; m < M2; m++) {
        mask2[n * M2 + m] ^= ( (m + n) % 3 == 0 );
      }
    }
  }
  // Expected
  // M1: eliminate based on Q. Or loop from q%3 to M1 step by 3
  // M2: total is ~1/3 N
  // Y sign fine
  for (unsigned int q = 0; q < Q; q++) {
    // 0, 1, 2 -> 0, 2, 1
    const size_t m1_start = (3 - q%3)%3;
    for (unsigned int d = 0; d < D; d++) {
      // x choice
      auto e = (d + q) % 2 ? expected1 : expected0;
      // M1, only when m + q = 0 mod 3
      for (unsigned int m1 = m1_start; m1 < M1; m1 += 3) {
        // M2 split into thirds
        for (unsigned int m2 = 0; m2 < M2; m2++) {
          const size_t idx = ((q * M1 + m1) * M2 + m2) * D + d;
          // 1/3 plus remainders
          e[idx] = ((3 * M2 - m2)%3 < (N%3)) + N/3;
          // y sign
          if (((d>>1) + q) % 2 == 1) e[idx] *= -1;
        }
      }
    }
  }

  if (use_ot_version) {
    store->check_DaBits(N * Q * D);
  } else {
    // Can also upper bound 4 * N Q M1 M2 D probably
    store->check_BoolTriples(N * Q * (3 * M1 * M2 * D + M1 * M2 + D));
    store->check_DaBits(4 * N * Q * M1 * M2 * D);
  }

  // Run
  auto start = clock_start();
  std::cout << "clock start" << std::endl;
  int64_t sent_bytes = 0;
  if (use_ot_version) {
    fmpz_t* y_p; new_fmpz_array(&y_p, N * Q * D);
    sent_bytes += store->b2a_single(N * Q * D, y, y_p);
    std::cout << "heavy mask: b2a timing " << sec_from(start) << std::endl;
    auto start2 = clock_start();
    bool* const combined_mask = new bool[N * Q * M1 * M2];
    sent_bytes += store->multiply_BoolShares_cross(
        N, Q * M1, M2, mask1, mask2, combined_mask);
    std::cout << "heavy mask: cross timing " << sec_from(start2) << std::endl;
    start2 = clock_start();
    sent_bytes += store->heavy_convert_mask(N, Q, M1 * M2, D,
        x, y_p, combined_mask, valid, bucket0, bucket1);
    std::cout << "heavy mask: convert timing " << sec_from(start2) << std::endl;
    delete[] combined_mask;
    clear_fmpz_array(y_p, N * Q * D);
  } else {
    // Expand mask2 across Q (copy for each Q).
    bool* const mask2_ext = new bool[N * Q * M2];
    for (unsigned int i = 0; i < N * Q; i++)
      memcpy(&mask2_ext[i * M2], &mask2[(i / Q) * M2], M2);
    sent_bytes += store->heavy_convert_mask_two(N, Q, M1, M2, D,
        x, y, mask1, mask2_ext, valid, bucket0, bucket1);
    delete[] mask2_ext;
    std::cout << "heavy mask: convert timing " << sec_from(start) << std::endl;
  }
  std::cout << "send_bytes: " << sent_bytes << std::endl;

  // Recombine / test
  if (server_num == 0) {
    fmpz_t* bucket0_other; new_fmpz_array(&bucket0_other, Q * M1 * M2 * D);
    fmpz_t* bucket1_other; new_fmpz_array(&bucket1_other, Q * M1 * M2 * D);
    recv_fmpz_batch(serverfd, bucket0_other, Q * M1 * M2 * D);
    recv_fmpz_batch(serverfd, bucket1_other, Q * M1 * M2 * D);
    fmpz_t tmp0; fmpz_init(tmp0);
    fmpz_t tmp1; fmpz_init(tmp1);
    for (unsigned int i = 0; i < Q * M1 * M2 * D; i++) {
      [[maybe_unused]] const size_t q = i / (M1 * M2 * D);
      [[maybe_unused]] const size_t m1 = (i / (M2 * D)) % M1;
      [[maybe_unused]] const size_t m2 = (i / D) % M2;
      [[maybe_unused]] const size_t d = i % D;

      fmpz_mod_add(tmp0, bucket0[i], bucket0_other[i], mod_ctx);
      fmpz_mod_add(tmp1, bucket1[i], bucket1_other[i], mod_ctx);

      // std::cout << "bucket[" << q << ", " << m1 << ", " << m2 << ", " << d << "]: ";
      // std::cout << "0 = " << fmpz_mod_get_signed(tmp0, Int_Modulus);
      // std::cout << ", expected = " << expected0[i] << ", \t";
      // std::cout << "1 = " << fmpz_mod_get_signed(tmp1, Int_Modulus);
      // std::cout << ", expected = " << expected1[i] << std::endl;

      assert(fmpz_mod_get_signed(tmp0, Int_Modulus) == expected0[i]);
      assert(fmpz_mod_get_signed(tmp1, Int_Modulus) == expected1[i]);
    }
    fmpz_clear(tmp0);
    fmpz_clear(tmp1);
    clear_fmpz_array(bucket0_other, Q * M1 * M2 * D);
    clear_fmpz_array(bucket1_other, Q * M1 * M2 * D);
  } else {
    send_fmpz_batch(serverfd, bucket0, Q * M1 * M2 * D);
    send_fmpz_batch(serverfd, bucket1, Q * M1 * M2 * D);
  }

  delete[] x;
  delete[] y;
  delete[] valid;
  delete[] mask1;
  delete[] mask2;
  clear_fmpz_array(bucket0, Q * M1 * M2 * D);
  clear_fmpz_array(bucket1, Q * M1 * M2 * D);
}

void test_abs_cmp(
    const int server_num, const int serverfd, CorrelatedStore* store) {
  // setup
  const size_t N = 10;
  fmpz_t* val0; new_fmpz_array(&val0, N);
  fmpz_t* val1; new_fmpz_array(&val1, N);

  // "Very small" compared to modulus
  const size_t range = 1000;
  fmpz_t range_f; fmpz_init_set_si(range_f, range);
  for (unsigned int i = 0; i < N; i++) {
    fmpz_randm(val0[i], seed, range_f);
    fmpz_mod_sub_si(val0[i], val0[i], range / 2, mod_ctx);

    fmpz_randm(val1[i], seed, range_f);
    fmpz_mod_sub_si(val1[i], val1[i], range / 2, mod_ctx);
  }
  // fmpz_set_si(val0[0], -(100 + 10*server_num));
  // fmpz_mod(val0[0], val0[0], Int_Modulus);
  // fmpz_set_si(val1[0], (10 + server_num));
  // fmpz_mod(val1[0], val1[0], Int_Modulus);

  // run
  fmpz_t* larger; new_fmpz_array(&larger, N);
  store->abs_cmp(N, val0, val1, larger);

  // test
  if (server_num == 0) {
    fmpz_t* larger_other; new_fmpz_array(&larger_other, N);
    recv_fmpz_batch(serverfd, larger_other, N);
    fmpz_t* val0_other; new_fmpz_array(&val0_other, N);
    fmpz_t* val1_other; new_fmpz_array(&val1_other, N);
    recv_fmpz_batch(serverfd, val0_other, N);
    recv_fmpz_batch(serverfd, val1_other, N);

    fmpz_t v0; fmpz_init(v0);
    fmpz_t v1; fmpz_init(v1);
    fmpz_t half; fmpz_init(half);
    fmpz_cdiv_q_ui(half, Int_Modulus, 2);

    fmpz_t v0_tmp; fmpz_init(v0_tmp);
    fmpz_t v1_tmp; fmpz_init(v1_tmp);

    for (unsigned int i = 0; i < N; i++) {
      fmpz_mod_add(larger[i], larger[i], larger_other[i], mod_ctx);
      // std::cout << "larger[" << i << "] = " << fmpz_get_ui(larger[i]) << std::endl;

      fmpz_mod_add(v0, val0[i], val0_other[i], mod_ctx);
      fmpz_set(v0_tmp, v0);

      fmpz_mod_add(v1, val1[i], val1_other[i], mod_ctx);
      fmpz_set(v1_tmp, v1);

      if (fmpz_cmp(v0_tmp, half) > 0) {  // > N/2, negate (abs)
        fmpz_sub(v0_tmp, Int_Modulus, v0);
      }
      if (fmpz_cmp(v1_tmp, half) > 0) {  // > N/2, negate (abs)
        fmpz_sub(v1_tmp, Int_Modulus, v1_tmp);
      }
      bool actual = fmpz_cmp(v0_tmp, v1_tmp) < 0;

      // std::cout << i << ": |" << fmpz_mod_get_signed(v0, Int_Modulus);
      // std::cout << (fmpz_is_one(larger[i]) ? "| < |" : "| > |");
      // std::cout << fmpz_mod_get_signed(v1, Int_Modulus);
      // std::cout << "|, \tactual: " << (actual ? "<" : ">") << std::endl;
      // std::cout << "(" << fmpz_get_ui(val0[i]) << " + " << fmpz_get_ui(val0_other[i]) << ") = ";
      // std::cout << fmpz_mod_get_signed(v0, Int_Modulus) << " vs ";
      // std::cout << fmpz_mod_get_signed(v1, Int_Modulus) << " = (";
      // std::cout << fmpz_get_ui(val1[i]) << " + " << fmpz_get_ui(val1_other[i]) << ")\n";
      assert(fmpz_is_one(larger[i]) == actual);
    }
    fmpz_clear(v0);
    fmpz_clear(v1);
    fmpz_clear(half);
    clear_fmpz_array(larger_other, N);
    clear_fmpz_array(val0_other, N);
    clear_fmpz_array(val1_other, N);
  } else {
    send_fmpz_batch(serverfd, larger, N);
    send_fmpz_batch(serverfd, val0, N);
    send_fmpz_batch(serverfd, val1, N);
  }

  clear_fmpz_array(larger, N);
  clear_fmpz_array(val0, N);
  clear_fmpz_array(val1, N);
}

void test_cmp_bit(
    const int server_num, const int serverfd, CorrelatedStore* store) {
  // Setup values
  const size_t bits = 2;
  const size_t max = 1ULL << bits;
  const size_t N = max * max;
  fmpz_t* x_bits; new_fmpz_array(&x_bits, N * bits);
  fmpz_t* y_bits; new_fmpz_array(&y_bits, N * bits);
  if (server_num == 0) {
    fmpz_t* x_bits_other; new_fmpz_array(&x_bits_other, N * bits);
    fmpz_t* y_bits_other; new_fmpz_array(&y_bits_other, N * bits);
    for (unsigned int x = 0; x < max; x++) {
      for (unsigned int y = 0; y < max; y++) {
        // std::cout << "case x = " << x << ", y = " << y << std::endl;
        // if (x == y) continue;  // should be 0, skip for smaller logging for now
        for (unsigned int b = 0; b < bits; b++) {
          size_t idx = bits * (x * max + y) + b;
          fmpz_randm(x_bits[idx], seed, Int_Modulus);
          fmpz_sub(x_bits_other[idx], Int_Modulus, x_bits[idx]);
          fmpz_add_ui(x_bits_other[idx], x_bits_other[idx], (x >> b) % 2);

          fmpz_randm(y_bits[idx], seed, Int_Modulus);
          fmpz_sub(y_bits_other[idx], Int_Modulus, y_bits[idx]);
          fmpz_add_ui(y_bits_other[idx], y_bits_other[idx], (y >> b) % 2);

          // std::cout << " bit " << b << ", idx = " << idx << std::endl;
          // std::cout << "  x bit: " << ((x >> b) % 2) << " = " << fmpz_get_ui(x_bits[idx]);
          // std::cout << " + " << fmpz_get_ui(x_bits_other[idx]) << std::endl;
          // std::cout << "  y bit: " << ((y >> b) % 2) << " = " << fmpz_get_ui(y_bits[idx]);
          // std::cout << " + " << fmpz_get_ui(y_bits_other[idx]) << std::endl;
        }
      }
    }
    send_fmpz_batch(serverfd, x_bits_other, N * bits);
    send_fmpz_batch(serverfd, y_bits_other, N * bits);
    clear_fmpz_array(x_bits_other, N * bits);
    clear_fmpz_array(y_bits_other, N * bits);
  } else {
    recv_fmpz_batch(serverfd, x_bits, N * bits);
    recv_fmpz_batch(serverfd, y_bits, N * bits);
  }

  store->check_Triples(3 * N * bits);

  // Eval cmp
  std::cout << "running" << std::endl;
  fmpz_t* ans; new_fmpz_array(&ans, N);
  store->cmp_bit(N, bits, x_bits, y_bits, ans);

  // Check values
  if (server_num == 0) {
    fmpz_t* ans_other; new_fmpz_array(&ans_other, N);
    recv_fmpz_batch(serverfd, ans_other, N);

    fmpz_t tmp; fmpz_init(tmp);
    for (unsigned int x = 0; x < max; x++) {
      for (unsigned int y = 0; y < max; y++) {
        size_t idx = x * max + y;
        fmpz_mod_add(tmp, ans[idx], ans_other[idx], mod_ctx);
        std::cout << x << (fmpz_is_one(tmp) ? " <  " : " >= " ) << y;
        std::cout << " (" << fmpz_get_ui(ans[idx]) << " + " << fmpz_get_ui(ans_other[idx]);
        std::cout << " = " << fmpz_get_ui(tmp) << ") " << std::endl;
        assert(fmpz_is_one(tmp) == (x < y));
      }
    }
    clear_fmpz_array(ans_other, N);
    fmpz_clear(tmp);
  } else {
    send_fmpz_batch(serverfd, ans, N);
  }

  clear_fmpz_array(ans, N);
  clear_fmpz_array(x_bits, N * bits);
  clear_fmpz_array(y_bits, N * bits);
}

void test_rand_bitshare(
    const int server_num, const int serverfd, CorrelatedStore* store) {
  const size_t N = 10;
  fmpz_t M; fmpz_init_set(M, Int_Modulus);
  const size_t b = fmpz_clog_ui(Int_Modulus, 2);

  // Expose get_BitShare
  class Test : CorrelatedStore {
  public:
    const BitShare* test_get() {return get_BitShare();};
  };

  store->check_BitShares(N);
  ((PrecomputeStore*)store)->print_Sizes();

  fmpz_t* buff; new_fmpz_array(&buff, N * (b + 1));
  for (unsigned int i = 0; i < N; i++) {
    const BitShare* x = ((Test*)store)->test_get();
    fmpz_set(buff[i * (b+1)], x->r);
    copy_fmpz_array(&buff[i * (b+1) + 1], x->rB, b);
    delete x;
  }
  reveal_fmpz_batch(serverfd, buff, N * (b + 1));

  fmpz_t val; fmpz_init(val);
  fmpz_t bit_val; fmpz_init(bit_val);

  // std::cout << "M = " << fmpz_get_ui(M) << " (" << b << " bits)" << std::endl;

  for (unsigned int i = 0; i < N; i++) {
    const size_t offset = i * (b + 1);
    fmpz_set(val, buff[offset]);
    // if (server_num == 0) std::cout << "val r = " << fmpz_get_ui(val) << std::endl;
    fmpz_zero(bit_val);
    for (int j = 0; j < b; j++) {
      // std::cout << " bit " << j << " = " << fmpz_get_ui(buff[offset + 1 + j]) << std::endl;
      fmpz_addmul_ui(bit_val, buff[offset + 1 + j], 1ULL << j);
    }
    // if (server_num == 0) std::cout << " bit val: " << fmpz_get_ui(bit_val) << std::endl;
    // Ensure rB represent r
    assert(fmpz_equal(val, bit_val));
    // Ensure r < M
    assert(fmpz_cmp(val, M) < 0);
  }
  clear_fmpz_array(buff, N  * (b + 1));
  fmpz_clear(val);
  fmpz_clear(bit_val);
  fmpz_clear(M);
}

void test_sign(
    const int server_num, const int serverfd, CorrelatedStore* store) {
  const size_t N = 10;

  fmpz_t* x; new_fmpz_array(&x, N);
  fmpz_t* actual; new_fmpz_array(&actual, N);
  fmpz_t* other; new_fmpz_array(&other, N);

  if (server_num == 0) {
    for (unsigned int i = 0; i < N; i++) {
      fmpz_randm(x[i], seed, Int_Modulus);
      fmpz_randm(other[i], seed, Int_Modulus);
      fmpz_mod_add(actual[i], x[i], other[i], mod_ctx);
      // std::cout << "val[" << i << "] = " << fmpz_get_ui(actual[i]) << " = ";
      // std::cout << fmpz_get_ui(x[i]) << " + " << fmpz_get_ui(other[i]) << std::endl;
    }
    send_fmpz_batch(serverfd, other, N);
  } else {
    recv_fmpz_batch(serverfd, x, N);
  }

  fmpz_t* s; new_fmpz_array(&s, N);
  // store->LSB(N, x, s);
  store->is_negative(N, x, s);

  if (server_num == 0) {
    recv_fmpz_batch(serverfd, other, N);

    fmpz_t tmp; fmpz_init(tmp);
    for (unsigned int i = 0; i < N; i++) {
      fmpz_mod_add(tmp, s[i], other[i], mod_ctx);
      // TODO: assert
      // std::cout << "val[" << i << "] = " << fmpz_get_ui(actual[i]);
      // std::cout << " has least sig bit " << fmpz_get_ui(tmp) << " = " << fmpz_get_ui(s[i]);
      // std::cout << " + " << fmpz_get_ui(other[i]) << std::endl;
      std::cout << "val[" << i << "] = " << fmpz_mod_get_signed(actual[i], Int_Modulus) << " (";
      std::cout << fmpz_get_ui(actual[i]) << ") is ";
      std::cout << (fmpz_is_one(tmp) == 1 ? "negative" : "non-negative") << std::endl;
    }

    fmpz_clear(tmp);
  } else {
    send_fmpz_batch(serverfd, s, N);
  }

  clear_fmpz_array(s, N);
  clear_fmpz_array(x, N);
  clear_fmpz_array(actual, N);
  clear_fmpz_array(other, N);
}

void runServerTest(const int server_num, const int serverfd) {

  // init
  OT_Wrapper* ot0 = new OT_Wrapper(server_num == 0 ? nullptr : "127.0.0.1", 60051);
  OT_Wrapper* ot1 = new OT_Wrapper(server_num == 1 ? nullptr : "127.0.0.1", 60052);
  PrecomputeStore* store = new PrecomputeStore(serverfd, server_num, ot0, ot1, batch_size, lazy);

  // store->maybe_Update();
  // std::cout << std::endl;

  // random adjusting. different numbers adjust seed.
  // Different per server to desync it
  fmpz_t tmp;
  fmpz_init(tmp);
  size_t rand_adjust;
  rand_adjust = (server_num == 0 ? 1 : 5);
  for (unsigned int i = 0; i < rand_adjust; i++)
    fmpz_randm(tmp, seed, Int_Modulus);
  fmpz_clear(tmp);

  test_HeavyConvert(server_num, serverfd, store);
  test_HeavyConvertMask_one(server_num, serverfd, store);
  test_HeavyConvertMask_two(server_num, serverfd, store);

  // Lazy 1 rand bitshare doesn't work right for some reason

  // test_cmp_bit(server_num, serverfd, store);
  test_rand_bitshare(server_num, serverfd, store);
  // test_sign(server_num, serverfd, store);
  test_abs_cmp(server_num, serverfd, store);

  delete ot0;
  delete ot1;
  delete store;
}

void serverTest() {
  std::cout << "Running server test" << std::endl;
  int sockfd = init_receiver();

  pid_t pid = fork();
  if (pid == 0) {
    int cli_sockfd = init_sender();

    runServerTest(0, cli_sockfd);
    close(cli_sockfd);
  } else {
    int newsockfd = accept_receiver(sockfd);

    // alter randomness to be different from the sender
    fmpz_t tmp;
    fmpz_init(tmp);
    size_t rand_adjust = 4;
    for (unsigned int i = 0; i < rand_adjust; i++)
      fmpz_randm(tmp, seed, Int_Modulus);
    fmpz_clear(tmp);

    runServerTest(1, newsockfd);

    close(newsockfd);
  }

  close(sockfd);
}

// 00, 01 = 0, 10 = 1, 11 = 0 (-1)
int main(int argc, char* argv[]) {
  init_constants();

  std::cout << "Using " << (use_ot_version ? "" : "non-") << "OT heavy convert\n";

  int server_num = -1;
  if (argc >= 2){
    server_num = atoi(argv[1]);
  }
  if (server_num == -1) {
    serverTest();
  } else if (server_num == 0) {
    int sockfd = init_receiver();
    int newsockfd = accept_receiver(sockfd);
    runServerTest(0, newsockfd);
    close(newsockfd);
    close(sockfd);
  } else if (server_num == 1) {
    int cli_sockfd = init_sender();
    runServerTest(1, cli_sockfd);
    close(cli_sockfd);
  }

  clear_constants();
}
