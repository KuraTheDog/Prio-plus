#include "share.h"

#include "constants.h"
#include "fmpz_utils.h"

extern "C" {
  #include "flint/flint.h"
  #include "flint/fmpz.h"
};

void init_client_packet(ClientPacket &p, const int N, const int NumMulInpGates) {
    p = new client_packet();

    p->N = N;
    p->NWires = NumMulInpGates;

    new_fmpz_array(&p->WireShares, NumMulInpGates);
    
    fmpz_init(p->f0_s);
    fmpz_init(p->g0_s);
    fmpz_init(p->h0_s);

    new_fmpz_array(&p->h_points, N);

    p->triple_share = new BeaverTripleShare();
}

void SplitShare(const fmpz_t val, fmpz_t A, fmpz_t B) {
    fmpz_randm(A,seed,Int_Modulus);
    fmpz_sub(B,val,A);
    fmpz_mod(B,B,Int_Modulus);
}

// Unused?
void SplitShare(const fmpz_t val, fmpz_t A, fmpz_t B, const int num_bits) {
    // num_bits < 32
    uint32_t mod = 1 << num_bits;

    fmpz_t max_val;
    fmpz_init(max_val);
    fmpz_set_ui(max_val, mod);
    fmpz_randm(A,seed,max_val);
    fmpz_xor(B,A,val);
    fmpz_clear(max_val);
}

BeaverTriple* NewBeaverTriple() {
    BeaverTriple* out = new BeaverTriple();

    fmpz_randm(out->A,seed,Int_Modulus);
    fmpz_randm(out->B,seed,Int_Modulus);
    fmpz_mul(out->C,out->A,out->B);
    fmpz_mod(out->C,out->C,Int_Modulus);

    return out;
}

BeaverTripleShare* BeaverTripleShares(const BeaverTriple* inp) {
    BeaverTripleShare* out = new BeaverTripleShare[2];

    SplitShare(inp->A,out[0].shareA,out[1].shareA);
    SplitShare(inp->B,out[0].shareB,out[1].shareB);
    SplitShare(inp->C,out[0].shareC,out[1].shareC);

    return out;
}
