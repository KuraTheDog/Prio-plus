typedef struct BitShare{
    char pk[32];
    unsigned int val;
    char signature[32];
} BitShare;

typedef struct IntShare{
    char pk[32];
    unsigned int val;
    char signature[32];
};

enum messageType{
    BIT_SUM,
    INIT_BIT_SUM,
    INT_SUM,
    INIT_INT_SUM
};

typedef struct initMsg{
    messageType type;
    int num_of_inputs;
}initMsg;


