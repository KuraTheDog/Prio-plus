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
    uint32_t* arr;
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
};

struct initMsg {
    messageType type;
    unsigned int num_of_inputs;
    unsigned int max_inp;
};

#endif
