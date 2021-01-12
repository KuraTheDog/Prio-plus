#ifndef TYPES_H
#define TYPES_H

#define PK_LENGTH 32

struct BitShare {
    char pk[PK_LENGTH];
    bool val;
    char signature[PK_LENGTH];
};

// For INT_SUM, AND_OP, OR_OP
struct IntShare {
    char pk[PK_LENGTH];
    uint32_t val;
    char signature[PK_LENGTH];
};

// For Max, Min
struct MaxShare {
    uint32_t* arr;
    char pk[PK_LENGTH];
    char signature[PK_LENGTH];
};

// For Var, Stddev
struct VarShare {
    char pk[PK_LENGTH];
    uint32_t val;
    uint32_t val_squared;
    char signature[PK_LENGTH];
};

struct LinRegShare {
    char pk[PK_LENGTH];
    uint32_t num_fields;
    uint32_t* vals;
    char signature[PK_LENGTH];
};

enum messageType {
    BIT_SUM,
    INT_SUM,
    AND_OP,
    OR_OP,
    MAX_OP,
    MIN_OP,
    VAR_OP,
    STDDEV_OP,
};

struct initMsg {
    messageType type;
    int num_of_inputs;
    int max_inp;
};

#endif
