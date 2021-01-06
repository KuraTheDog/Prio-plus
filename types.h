#ifndef TYPES_H
#define TYPES_H

#define PK_LENGTH 32

struct BitShare{
    char pk[PK_LENGTH];
    unsigned int val;
    char signature[PK_LENGTH];
};

struct IntShare{
    char pk[PK_LENGTH];
    unsigned int val;
    char signature[PK_LENGTH];
};

struct AndShare{
    char pk[PK_LENGTH];
    unsigned int val;
    char signature[PK_LENGTH];
};

struct OrShare{
    char pk[PK_LENGTH];
    unsigned int val;
    char signature[PK_LENGTH];
};

struct MaxShare{
    uint32_t* arr;
    char pk[PK_LENGTH];
    char signature[PK_LENGTH];
};

struct VarShare {
    char pk[PK_LENGTH];
    unsigned int val;
    unsigned int val_squared;
    char signature[PK_LENGTH];
};

enum messageType{
    BIT_SUM,
    INT_SUM,
    AND_OP,
    OR_OP,
    MAX_OP,
    VAR_OP,
};

struct initMsg{
    messageType type;
    int num_of_inputs;
    int max_inp;
};

#endif
