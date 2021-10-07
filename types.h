#ifndef TYPES_H
#define TYPES_H

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
    COUNTMIN_OP,
    HEAVY_OP,
};

struct initMsg {
    messageType type;
    unsigned int num_bits;
    unsigned int num_of_inputs;
    unsigned int max_inp;
};

struct HeavyConfig {
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

#endif
