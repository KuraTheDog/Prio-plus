#ifndef TYPES_H
#define TYPES_H

#define TAG_LENGTH 32

struct BitShare {
    char tag[TAG_LENGTH];
    bool val;
};

// For INT_SUM, AND_OP, OR_OP
struct IntShare {
    char tag[TAG_LENGTH];
    uint64_t val;
};

// For Max, Min
struct MaxShare {
    char tag[TAG_LENGTH];
    uint64_t* arr;
};

// For Var, Stddev
struct VarShare {
    char tag[TAG_LENGTH];
    uint64_t val;
    uint64_t val_squared;
};

struct LinRegShare {
    char tag[TAG_LENGTH];
    uint64_t* x_vals;   // Feats: x0, x1, ... d-1
    uint64_t y;         // Target: y 1
    uint64_t* x2_vals;  // Quadratic in feats: x0^2, x0 x1, x1^2, ... d(d-1)/2
    uint64_t* xy_vals;  // Feat * target: x0 y, x1 y, (d-1)
};

struct FreqShare {
    char tag[TAG_LENGTH];
    bool* arr;  // Could theoretically be compacted
};

// Multiple frequency arrays of varying lengths.
struct MultiFreqShare {
    char tag[TAG_LENGTH];
    bool** arr;
};

enum messageType {
    NONE_OP,
    BIT_SUM_OP,
    INT_SUM_OP,
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

#endif
