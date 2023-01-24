#ifndef TYPES_H
#define TYPES_H

#include <iostream>
#include "utils.h"

#define PK_LENGTH 32

struct BitShare {
    char pk[PK_LENGTH];
    bool val;
};

// For INT_SUM, AND_OP, OR_OP
struct IntShare {
    char pk[PK_LENGTH];
    uint64_t val;
};

// For Max, Min
struct MaxShare {
    char pk[PK_LENGTH];
    uint64_t* arr;
};

// For Var, Stddev
struct VarShare {
    char pk[PK_LENGTH];
    uint64_t val;
    uint64_t val_squared;
};

struct LinRegShare {
    char pk[PK_LENGTH];
    uint64_t* x_vals;   // Feats: x0, x1, ... d-1
    uint64_t y;         // Target: y 1
    uint64_t* x2_vals;  // Quadratic in feats: x0^2, x0 x1, x1^2, ... d(d-1)/2
    uint64_t* xy_vals;  // Feat * target: x0 y, x1 y, (d-1)
};

struct FreqShare {
    char pk[PK_LENGTH];
    bool* arr;  // Could theoretically be compacted
};

// Multiple frequency arrays of varying lengths.
struct MultiFreqShare {
    char pk[PK_LENGTH];
    bool** arr;
};

enum messageType {
    NONE_OP,
    BIT_SUM,
    INT_SUM,
    AND_OP,
    OR_OP,
    MAX_OP,
    MIN_OP,
    VAR_OP,
    STDDEV_OP,
    LINREG_OP,
    FREQ_OP,
    HEAVY_OP,
    MULTI_HEAVY_OP,
};

struct initMsg {
    messageType type;
    unsigned int num_bits;
    unsigned int num_of_inputs;
    unsigned int max_inp;
};


struct CountMinConfig {
    /* Target parameters
    Return all heavy > t
    With odds (1 - delta), nothing < (1 - eps) t is returned
    */
    double t;  // Target heavy. 0.1 is >10%
    // double eps;
    // double delta;

    /* Tree params
    For full tree eps/delta, per layer eps/delta values
    Since recersive, can be slightly more efficient
    */
    // double eps2;  // 1-eps = (1 - eps)^L, so eps' = 1 - log_L (1 - eps)
    // double delta2;  // delta' = delta (1 - eps) t / 2 L
    size_t L;  // # layers, input_bits - ceil(log_2(w d))

    /* Input parameters
    w = a / eps t, d = log_a(1 / delta)
    For tree, use eps2, delta2 in place
    */
    size_t w;  // hash range
    size_t d;  // # hashes
};

struct MultiHeavyConfig {
    // n = 2^num_bits, m = numreqs
    // Threshold K. find top K
    // Delta: Failure (0.05)
    // Q = log_2(1/delta) copies of single-heavy
    // C: Return among top CK. small > 2 + (2/K ln(1/delta))^(1/3), e.g. 10
    // C' >= 1, < (C-(2/K ln(1/delta))^(1/3)) / 2 : min random subset size
    // R = layers. "log n", but I think 1 works.
    // B = hash range, (CK)^2 / (2 ln 2)
    // Also should have B < input size
    // Countmin parameters

    // root values. The rest can be derived from these
    const size_t K;
    const double delta;
    const unsigned int delta_inv; // = (unsigned int) floor(1/delta);
    const size_t Q; // = LOG2(delta_inv);
    const size_t C_prime = 1;
    const size_t C = 2 * C_prime + 1;  // +1 fine for delta >= 0.00034
    const size_t R = 1;  // Currently not implemented
    const double twoln2 = 1.3862943611198906188;
    const size_t B; // = ceil((C * K) * (C * K) / twoln2);

    // Q parallel iterations.
    // Each splits input into B substreams.
    // Each substream runs a SH instance. For anti-parallel, this is 2 log n
    // QB total SH instances.

    const size_t num_bits;
    const size_t SH_depth;

    // TODO: Countmin
    CountMinConfig count_min_config;

    MultiHeavyConfig(size_t K, double delta, size_t num_bits)
    : K(K)
    , delta(delta)
    , delta_inv((unsigned int) floor(1/delta))
    , Q(LOG2(delta_inv))
    , B(ceil((C * K) * (C * K) / twoln2))
    , num_bits(num_bits)
    , SH_depth(2 * num_bits)
    {}

    void print() {
        std::cout << "Multi Heavy Params: " << std::endl;
        std::cout << "\tK = # heavy target: " << K << std::endl;
        std::cout << "\tdelta = max failure: " << delta << std::endl;
        std::cout << "\tQ = # SH : " << Q << std::endl;
        std::cout << "\tC = among top CK: " << C << std::endl;
        std::cout << "\tR = # SH layers: " << R << std::endl;
        std::cout << "\tB = substream hash range: " << B << std::endl;
        std::cout << "\tnum_bits: " << num_bits << std::endl;
        std::cout << "\tSH depth: " << SH_depth << std::endl;
        // TODO: also count min
    }
};

#endif
