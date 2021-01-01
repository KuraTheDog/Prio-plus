#ifndef TYPES_H
#define TYPES_H

// TODO: 32 as global constant.

struct BitShare{
    char pk[32];
    unsigned int val;
    char signature[32];
};

struct IntShare{
    char pk[32];
    unsigned int val;
    char signature[32];
};

struct AndShare{
    char pk[32];
    unsigned int val;
    char signature[32];
};

struct OrShare{
    char pk[32];
    unsigned int val;
    char signature[32];
};

struct MaxShare{
    uint32_t* arr;
    char pk[32];
    char signature[32];
};

struct VarShare {
    char pk[32];
    unsigned int val;
    unsigned int val_squared;
    // TODO: snips
    char signature[32];
};

enum messageType{
    BIT_SUM,
    INIT_BIT_SUM,
    INT_SUM,
    INIT_INT_SUM,
    AND_OP,
    INIT_AND_OP,
    OR_OP,
    INIT_OR_OP,
    MAX_OP,
    INIT_MAX_OP,
    VAR_OP,
    INIT_VAR_OP
};

struct initMsg{
    messageType type;
    int num_of_inputs;
    int max_inp;
};

#endif
