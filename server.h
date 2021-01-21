#ifndef SERVER_H
#define SERVER_H

#include "circuit.h"
#include "constants.h"
#include "fmpz_utils.h"
#include "net_share.h"
#include "share.h"

extern "C" {
    #include "poly/poly_batch.h"
    #include "poly/poly_once.h"
}

// Return type of different ops
enum returnType {
    RET_INVALID,    // Too many inputs invalid
    RET_ANS,        // Success, Returning ans. For the one thread that does the computation
    RET_NO_ANS,     // Success, no ans. For forking and support server.
};

// Only used in unused BatchPre.Interp
struct BatchPoly {
    fmpz_mod_poly_t fpoly;

    // Below 2 are unused in test_circuit.

    fmpz_t* Eval(const fmpz_t* xPointsIn, const int n) {
        return poly_batch_evaluate(fpoly, n, xPointsIn);
    }

    void EvalOnce(const fmpz_t x, fmpz_t out) {
        fmpz_init(out);
        poly_batch_evaluate_once(fpoly, x, out);
    }

    ~BatchPoly() {
        poly_batch_clear(fpoly);
    }
};

struct BatchPre {
    precomp_t pre;

    // Newbatch
    BatchPre(const fmpz_t *xPointsIn, const int n) {
        // std::cout << " new BatchPre , n = " << n << ", xPointsIn = [";
        // for (int i =0; i < n; i++) {
        //     if (i > 0) std::cout << ", ";
        //     fmpz_print(xPointsIn[i]);
        // }
        // std::cout << "]\n";
        poly_batch_precomp_init(&pre, Int_Modulus, n, xPointsIn);
    }

    // Unused.
    BatchPoly* Interp(const fmpz_t* yPointsIn, const int n) {
        BatchPoly* bpoly = new BatchPoly();

        poly_batch_init(bpoly->fpoly, &pre);

        poly_batch_interpolate(bpoly->fpoly, &pre, yPointsIn);

        return bpoly;
    }

    ~BatchPre() {
        poly_batch_precomp_clear(&pre);
    }
};

struct PreX {
    BatchPre *batchPre;
    precomp_x_t pre;

    // Replaces NewEvalPoint
    PreX(BatchPre* b, const fmpz_t x) {
        batchPre = b;
        precomp_x_init(&pre, &batchPre->pre, x);
    }

    ~PreX() {
        precomp_x_clear(&pre);
    }

    void Eval(const fmpz_t* yValues, fmpz_t out) {
        precomp_x_eval(&pre, yValues, out);
    }
};

struct CheckerPreComp {
    fmpz_t x;

    BatchPre *degN;
    BatchPre *deg2N;

    PreX *xN;
    PreX *x2N;

    CheckerPreComp(const Circuit* ckt) {
        int n = ckt->NumMulGates() + 1;
        int N = NextPowerofTwo(n);

        degN = new BatchPre(roots, N);
        deg2N = new BatchPre(roots2, 2 * N);
    }

    void setCheckerPrecomp(const fmpz_t val) {
        fmpz_init_set(x, val);

        xN = new PreX(degN, x);
        x2N = new PreX(deg2N, x);
    }
};

struct Checker {
    int serveridx;     // id of this server
    ClientPacket req;  // Client packet
    Circuit *ckt;      // Validation circuit

    int n;  // number of mult gates
    int N;  // NextPowerofTwo(n)

    fmpz_t *pointsF;  // Points on f. f(i) = ith mul gate left input
    fmpz_t *pointsG;  // Points on g. g(i) = ith mul gate right input
    fmpz_t *pointsH;  // points on h. Want to check if h = f * g

    // For sigma = [r * (f(r) * g(r) - h(r))]
    fmpz_t evalF;  // [f(r)]
    fmpz_t evalG;  // [r * g(r)]
    fmpz_t evalH;  // [r * h(r)]

    Checker(Circuit* c, const int idx) {
        ckt = c;
        serveridx = idx;
        n = ckt->NumMulGates()+1;
        N = NextPowerofTwo(n);

        new_fmpz_array(&pointsF, N);
        new_fmpz_array(&pointsG, N);
        new_fmpz_array(&pointsH, 2 * N);

        fmpz_init(evalF);
        fmpz_init(evalG);
        fmpz_init(evalH);
    }

    ~Checker() {
        clear_fmpz_array(pointsF, N);
        clear_fmpz_array(pointsG, N);
        clear_fmpz_array(pointsH, 2 * N);
        fmpz_clear(evalF);
        fmpz_clear(evalG);
        fmpz_clear(evalH);
    }

    void setReq(const ClientPacket pkt) {
        req = pkt;
        ckt->ImportWires(pkt, serveridx);
    }

    void evalPoly(const CheckerPreComp *pre) {
        // std::cout << "evalPoly" << std::endl;
        std::vector<Gate*> mulgates = ckt->MulGates();
        // Get constant terms from packet
        fmpz_set(pointsF[0], req->f0_s);
        fmpz_set(pointsG[0], req->g0_s);
        fmpz_set(pointsH[0], req->h0_s);

        // For all multiplication triples a_i * b_i = c_i
        //    polynomial [f(x)] has [f(i)] = [a_i]
        //    polynomial [g(x)] has [g(i)] = [b_i]
        for (int i = 1; i < n; i++) {
            fmpz_set(pointsF[i], mulgates[i-1]->ParentL->WireValue);
            fmpz_set(pointsG[i], mulgates[i-1]->ParentR->WireValue);
            // Set even values of h to be output wires.
            fmpz_set(pointsH[2*i], mulgates[i-1]->WireValue);
        }

        // Grab odd values of h from the packet.
        for (int j = 0; j < N; j++) {
            fmpz_set(pointsH[2 * j + 1], req->h_points[j]);
        }

        // set evals
        pre->xN->Eval(pointsF, evalF);
        pre->xN->Eval(pointsG, evalG);
        fmpz_mul(evalG, evalG, pre->x);
        fmpz_mod(evalG, evalG, Int_Modulus);
        pre->x2N->Eval(pointsH, evalH);
        fmpz_mul(evalH, evalH, pre->x);
        fmpz_mod(evalH, evalH, Int_Modulus);
    }

    CorShare* CorShareFn(const CheckerPreComp *pre) {
        evalPoly(pre);
        // std::cout << "CorShareFn" << std::endl;
        CorShare* out = new CorShare();

        fmpz_sub(out->shareD, evalF, req->triple_share->shareA);
        fmpz_mod(out->shareD, out->shareD, Int_Modulus);

        fmpz_sub(out->shareE, evalG, req->triple_share->shareB);
        fmpz_mod(out->shareE, out->shareE, Int_Modulus);

        return out;
    }

    Cor* CorFn(const CorShare* cs0, const CorShare* cs1) {
        Cor* out = new Cor();
        // std::cout << "CorFn" << std::endl;
        fmpz_add(out->D, cs0->shareD, cs1->shareD);
        fmpz_mod(out->D, out->D, Int_Modulus);

        fmpz_add(out->E, cs0->shareE, cs1->shareE);
        fmpz_mod(out->E, out->E, Int_Modulus);

        return out;
    }

    // To be fixed. Both servers need to use same random seed. Using constant 1 instead now.
    void randSum(fmpz_t out, const fmpz_t* arr) {
        int len = ckt->result_zero.size() + 1;

        fmpz_t tmp;
        fmpz_init(tmp);

        for (int i = 0; i < len; i++) {
            // fmpz_randm(tmp, seed, Int_Modulus);
            fmpz_set_ui(tmp, 1);
            fmpz_mul(tmp, tmp, arr[i]);
            fmpz_mod(tmp, tmp, Int_Modulus);

            fmpz_add(out, out, tmp);
        }

        fmpz_mod(out, out, Int_Modulus);

        fmpz_clear(tmp);
    }

    void OutShare(fmpz_t out, const Cor* corIn) {
        fmpz_t mulCheck;
        fmpz_t term;

        fmpz_init(mulCheck);
        fmpz_init(term);

        if (serveridx == 0) {
            fmpz_mul(mulCheck, corIn->D, corIn->E);
            fmpz_mod(mulCheck, mulCheck, Int_Modulus);
        }

        fmpz_mul(term, corIn->D, req->triple_share->shareB);
        fmpz_mod(term, term, Int_Modulus);
        fmpz_add(mulCheck, mulCheck, term);
        fmpz_mod(mulCheck, mulCheck, Int_Modulus);

        fmpz_mul(term, corIn->E, req->triple_share->shareA);
        fmpz_mod(term, term, Int_Modulus);
        fmpz_add(mulCheck, mulCheck, term);
        fmpz_mod(mulCheck, mulCheck, Int_Modulus);

        fmpz_add(mulCheck, mulCheck, req->triple_share->shareC);
        fmpz_mod(mulCheck, mulCheck, Int_Modulus);

        fmpz_sub(mulCheck, mulCheck, evalH);
        fmpz_mod(mulCheck, mulCheck, Int_Modulus);

        fmpz_t* arr;
        int num_zero_gates = ckt->result_zero.size();
        new_fmpz_array(&arr, num_zero_gates+1);

        fmpz_set(arr[0], mulCheck);

        for (int i = 0; i < num_zero_gates; i++) {
            fmpz_set(arr[i+1], ckt->result_zero[i]->WireValue);
        }

        randSum(out, arr);

        clear_fmpz_array(arr, num_zero_gates);
        fmpz_clear(mulCheck);
        fmpz_clear(term);
    }

    bool OutputIsValid(const fmpz_t output0, const fmpz_t output1) {
        fmpz_t out;
        fmpz_init(out);

        fmpz_add(out, output0, output1);
        fmpz_mod(out, out, Int_Modulus);

        bool ans = fmpz_is_zero(out);
        fmpz_clear(out);

        return ans;
    }

};

bool multiplyBoolShares(const int serverfd, const int server_num, const bool x, const bool y, const BooleanBeaverTriple triple) {
    bool z, d, e, d_this, e_this, d_other, e_other;

    d_this = x ^ triple.a;
    e_this = y ^ triple.b;

    send_bool(serverfd, d_this);
    recv_bool(serverfd, d_other);
    send_bool(serverfd, e_this);
    recv_bool(serverfd, e_other);

    d = d_this ^ d_other;
    e = e_this ^ e_other;

    z = triple.c ^ (x and e) ^ (y and d);
    if (server_num == 0)
        z ^= (d and e);

    return z;
}

void multiplyArithmeticShares(const int serverfd, const int server_num, const fmpz_t x, const fmpz_t y, fmpz_t z, const BeaverTriple* triple) {
    fmpz_t d, d_other, e, e_other;

    fmpz_init(d);
    fmpz_init(e);
    fmpz_init(d_other);
    fmpz_init(e_other);

    fmpz_add(d, x, triple->A);  // [x] - [a]
    fmpz_mod(d, d, Int_Modulus);
    fmpz_add(e, y, triple->B);  // [y] - [b]
    fmpz_mod(e, e, Int_Modulus);

    send_fmpz(serverfd, d);
    recv_fmpz(serverfd, d_other);
    send_fmpz(serverfd, e);
    recv_fmpz(serverfd, e_other);

    fmpz_add(d, d, d_other);  // x - a
    fmpz_mod(d, d, Int_Modulus);
    fmpz_add(e, e, e_other);  // y - b
    fmpz_mod(e, e, Int_Modulus);

    // [xy] = [c] + [x] e + [y] d - de
    // optimization: mod more often? mod between mul and add?
    fmpz_set(z, triple->C);
    fmpz_addmul(z, x, e);
    fmpz_addmul(z, y, d);
    if (server_num == 0)
        fmpz_submul(z, d, e);
    fmpz_mod(z, z, Int_Modulus);

    fmpz_clear(d);
    fmpz_clear(d_other);
    fmpz_clear(e);
    fmpz_clear(e_other);
}

// Convert bit share [x]_2 into [x]_p using a daBit.
void b2a_daBit(const int serverfd, const int server_num, const DaBit* dabit, const bool x, fmpz_t &xp) {

    const bool v_this = x ^ dabit->b2;
    bool v_other;
    send_bool(serverfd, v_this);
    recv_bool(serverfd, v_other);
    const bool v = v_this ^ v_other;

    // [x]_p = v + [b]_p - 2 v [b]_p. v only added for one server.
    // So we only add v on server 1 (when v = 1)
    if (v) {  // If v = 1, then [x]_p = (0/1) - [b]_p
        fmpz_set_ui(xp, server_num);
        fmpz_sub(xp, xp, dabit->bp);
        fmpz_mod(xp, xp, Int_Modulus);
    } else {  // If v = 0, then [x]_p = [b]_p
        fmpz_set(xp, dabit->bp);
    }
}

DaBit* generateDaBit(const int serverfd, const int server_num, const BeaverTriple* triple) {
    // Answer
    DaBit* dabit = new DaBit();

    // Create local
    DaBit* dabit0 = new DaBit();
    DaBit* dabit1 = new DaBit();
    makeLocalDaBit(dabit0, dabit1);

    // Exchange
    DaBit* tmp_dabit = new DaBit();
    if (server_num == 0) {
        send_DaBit(serverfd, dabit1);
        recv_DaBit(serverfd, tmp_dabit);
        delete dabit1;
        dabit1 = tmp_dabit;
    } else {
        send_DaBit(serverfd, dabit0);
        recv_DaBit(serverfd, tmp_dabit);
        delete dabit0;
        dabit0 = tmp_dabit;
    }

    // Xor boolean shares
    dabit->b2 = dabit0->b2 ^ dabit1->b2;

    // Xor arithmetic shares, using a xor b = a + b - 2ab
    fmpz_t z;
    fmpz_init(z);
    multiplyArithmeticShares(serverfd, server_num, dabit0->bp, dabit1->bp, z, triple);

    fmpz_add(dabit->bp, dabit0->bp, dabit1->bp);
    fmpz_submul_ui(dabit->bp, z, 2);
    fmpz_mod(dabit->bp, dabit->bp, Int_Modulus);

    fmpz_clear(z);
    delete dabit0;
    delete dabit1;

    return dabit;
}

EdaBit* generateEdaBit(const int serverfd, const int server_num, const size_t n, const BooleanBeaverTriple* triples, const DaBit* dabit){
    // Answer bit
    EdaBit* edabit = new EdaBit(n);

    // Create local
    EdaBit* edabit0 = new EdaBit(n);
    EdaBit* edabit1 = new EdaBit(n);
    makeLocalEdaBit(edabit0, edabit1, n);

    // Exchange
    EdaBit* tmp_edabit = new EdaBit(n);
    if (server_num == 0) {
        send_EdaBit(serverfd, edabit1, n);
        recv_EdaBit(serverfd, tmp_edabit, n);
        delete edabit1;
        edabit1 = tmp_edabit;
    } else {
        send_EdaBit(serverfd, edabit0, n);
        recv_EdaBit(serverfd, tmp_edabit, n);
        delete edabit0;
        edabit0 = tmp_edabit;
    }
    // Add arithmetic shares
    fmpz_add(edabit->r, edabit0->r, edabit1->r);
    fmpz_mod(edabit->r, edabit->r, Int_Modulus);

    // Add binary shares via circuit
    // c_{i+1} = c_i xor ((x_i xor c_i) and (y_i xor c_i))
    // output z_i = x_i xor y_i xor c_i
    bool carry = false;
    for (int i = 0; i < n; i++) {
        edabit->b[i] = carry ^ edabit0->b[i] ^ edabit1->b[i];
        bool x = carry ^ edabit0->b[i];
        bool y = carry ^ edabit1->b[i];
        bool z = multiplyBoolShares(serverfd, server_num, x, y, triples[i]);
        carry ^= z;
    }

    // Convert carry to arithmetic [b_n]_p
    fmpz_t bpn;
    fmpz_init(bpn);
    b2a_daBit(serverfd, server_num, dabit, carry, bpn);

    // Subtract out 2^n * [b_n]_p from r
    fmpz_submul_ui(edabit->r, bpn, 1 << n);
    fmpz_mod(edabit->r, edabit->r, Int_Modulus);

    fmpz_clear(bpn);
    delete edabit0;
    delete edabit1;

    return edabit;
}

#endif
