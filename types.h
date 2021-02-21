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
    uint32_t* arr;
    char pk[PK_LENGTH];
};

// For Var, Stddev
struct VarShare {
    char pk[PK_LENGTH];
    uint64_t val;
    uint64_t val_squared;
};

struct LinRegShare {
    char pk[PK_LENGTH];
    uint64_t* x_vals;  // x, x^2, x^3, ...
    uint64_t* y_vals;  // y, xy, x^2 y, ...
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
