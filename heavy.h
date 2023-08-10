#ifndef HEAVY_H
#define HEAVY_H

#include <iostream>
#include <math.h>  // For log

#include "fmpz_utils.h"
#include "hash.h"
#include "utils.h"

struct CountMinConfig {
  // d hashes from input bits to range w
  // Standard: w = ⌈e/ε⌉ and d = ⌈ln 1/δ⌉
  // Currently just assuming eps = delta

  // Always >= actual
  // less than eps * total over actual with probability at least 1 - delta

  const double delta; // odds outside of error range
  const double eps;   // additive max error
  const size_t d;     // num hashes
  const size_t w;     // Hash range

  CountMinConfig(const double delta, const double eps)
  : delta(delta)
  , eps(eps)
  , d(ceil(log(1. / delta)))
  , w(ceil(exp(1) / eps))
  {}

  void print() const {
    std::cout << "Count-Min Params: " << std::endl;
    std::cout << "\tdelta = max failure: " << delta << std::endl;
    std::cout << "\teps = additive max error: " << eps << std::endl;
    std::cout << "\td = num hashes: " << d << std::endl;
    std::cout << "\tw = Hash range : " << w << std::endl;
  }
};

struct MultiHeavyConfig {
  // n = 2^num_bits, m = numreqs
  // Threshold K. find top K
  // Delta: Failure (0.05)
  // Q = log_2(1/delta) copies of single-heavy
  // B = hash range, (C)^2 / (2 ln 2)
  // R = halving layers, log(n). R = 1 is just the always survive layer #0
  // Also should have B < input size, else more efficient to just do freq
  // Countmin parameters
  // TODO: Update with latest parameter computations

  // root values. The rest can be derived from these
  const size_t K;
  const double delta;
  const double eps;
  const unsigned int delta_inv;  // floor(1/delta), for maths
  const size_t Q;
  const double ln2_inv = 1.44269504;  // 1 / (ln 2)
  const size_t B;
  const size_t R;

  // Q parallel iterations.
  // Each splits input into B substreams.
  // Each substream runs a SH instance. For anti-parallel, this is 2 log n
  // QB total SH instances.

  const size_t num_bits;
  const size_t D;  // SH Depth

  const CountMinConfig countmin_cfg;

  MultiHeavyConfig(size_t K, double delta, size_t num_bits, double eps, size_t layers = 0)
  : K(K)
  , delta(delta)
  , eps(eps)
  , delta_inv((unsigned int) floor(1/delta))
  , Q(LOG2(delta_inv))
  , B(ceil(K * ln2_inv))
  , R(layers ? layers : num_bits)
  , num_bits(num_bits)
  , D(num_bits /*+ 1*/)  // +1 if extra constant term in SH for sanity check ec.
  , countmin_cfg(CountMinConfig(delta, eps))
  {}

  void print() const {
      std::cout << "Multi Heavy Params: " << std::endl;
      std::cout << "\tK = # heavy target: " << K << std::endl;
      std::cout << "\tdelta = max failure: " << delta << std::endl;
      std::cout << "\teps = tolerance: " << eps << std::endl;
      std::cout << "\tQ = # SH : " << Q << std::endl;
      std::cout << "\tB = substream hash range: " << B << std::endl;
      std::cout << "\tR = # halving layers: " << R << std::endl;
      std::cout << "\tnum_bits: " << num_bits << std::endl;
      std::cout << "\tD = SH depth: " << D << std::endl;
      countmin_cfg.print();
  }
};

/*
d x w counts struct.
d hashes range [w]
Add(x) -> counts[i, hash_i(x)] += 1 for each i in d
query(x) -> min_i[i, hash_i(x)]

Always >= actual
less than eps * total over actual with probability at least 1 - delta

Main use is aggregating via shares, so counts is exposed to be copied into.
Add is also generally unused, but kept for testing

TODO: Check w * d < input range.
Otherwise it's more space efficient to just to a normal frequency vector
Can either warn or actually replace with frequency.
*/

struct CountMin {
  const CountMinConfig cfg;

  // Outwardly accessible, so can be set (e.g. secret shared first)
  fmpz_t* counts = nullptr;    // Should be w * d
  HashStore* store = nullptr;  // Should be d hashes w output

  CountMin(const CountMinConfig cfg)
  : cfg(cfg)
  {}

  // Seperate from constructor, so can be copied into instead
  void setStore(const size_t input_bits, flint_rand_t hash_seed_arg) {
    store = new HashStorePoly(cfg.d, input_bits, cfg.w, hash_seed_arg, 2);
  }

  // Seperate from constructor, so can be copied into instead
  void init() {
    new_fmpz_array(&counts, cfg.w * cfg.d);
  }

  ~CountMin() {
    delete store;
    if (counts)
      clear_fmpz_array(counts, cfg.w * cfg.d);
  }

  void add(const unsigned int x, const unsigned int amount = 1) {
    if (!store) {
      std::cout << "ERROR: call countmin setStore first" << std::endl;
      return;
    }
    if (!counts) {
      std::cout << "ERROR: countmin set or init counts first" << std::endl;
      return;
    }
    fmpz_t out; fmpz_init(out);
    for (unsigned int i = 0; i < cfg.d; i++) {
      store->eval(i, x, out);
      unsigned int hx = fmpz_get_ui(out);
      // std::cout << "hash_" << i << "(" << x << ") = " << hx << std::endl;
      size_t idx = i * cfg.w + hx;
      fmpz_add_ui(counts[idx], counts[idx], amount);
      // std::cout << "adding " << amount << " to " << idx << " = [" << i << ", " << hx << "]" << std::endl;
    }
    fmpz_clear(out);
  }

  int query(const unsigned int x) {
    if (!store) {
      std::cout << "ERROR: call countmin setStore first" << std::endl;
      return -1;
    }
    if (!counts) {
      std::cout << "ERROR: countmin set or init counts first" << std::endl;
      return -1;
    }

    unsigned int min_count = UINT_MAX;

    fmpz_t out; fmpz_init(out);

    for (unsigned int i = 0; i < cfg.d; i++) {
      store->eval(i, x, out);
      unsigned int hx = fmpz_get_ui(out);
      unsigned int freq = fmpz_get_ui(counts[i * cfg.w + hx]);
      // std::cout << "h_" << i << "(" << x << ") = " << hx << ", freq = " << freq << std::endl;
      if (freq < min_count) {
        min_count = freq;
      }
    }

    fmpz_clear(out);

    return min_count;
  }

  void print() {
    for (unsigned int i = 0; i < cfg.d; i++) {
      for (unsigned int j = 0; j < cfg.w; j++) {
        uint64_t x = fmpz_get_ui(counts[i * cfg.w + j]);
        if (x)
          std::cout << x << " ";
        else
          std::cout << "_ ";
      }
      std::cout << std::endl;
    }
  }
};

#endif
