#ifndef SERVER_H
#define SERVER_H

#include "circuit.h"
#include "constants.h"
#include "fmpz_utils.h"

extern "C" {
    #include "poly/poly_batch.h"
    #include "poly/poly_once.h"
}

// Return type of different ops
enum returnType {
    RET_INVALID,    // Too many inputs invalid
    RET_ANS,        // Success, Returning ans. For server returning an answer
    RET_NO_ANS,     // Success, no ans. For support server.
};

// Only used in unused BatchPre.Interp
/*
struct BatchPoly {
    fmpz_mod_poly_t fpoly;

    // Below 2 are unused in test_circuit.

    fmpz_t* Eval(const fmpz_t* const xPointsIn, const int n) {
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
*/

struct BatchPre {
    precomp_t pre;

    // Newbatch
    BatchPre(const fmpz_t* const xPointsIn, const int n) {
        // std::cout << " new BatchPre , n = " << n << ", xPointsIn = [";
        // for (int i =0; i < n; i++) {
        //     if (i > 0) std::cout << ", ";
        //     fmpz_print(xPointsIn[i]);
        // }
        // std::cout << "]\n";
        poly_batch_precomp_init(&pre, Int_Modulus, n, xPointsIn);
    }

    // Unused.
    /*
    BatchPoly* Interp(const fmpz_t* const yPointsIn, const int n) {
        BatchPoly* bpoly = new BatchPoly();

        poly_batch_init(bpoly->fpoly, &pre);

        poly_batch_interpolate(bpoly->fpoly, &pre, yPointsIn);

        return bpoly;
    }
    */

    ~BatchPre() {
        poly_batch_precomp_clear(&pre);
    }
};

struct PreX {
    const BatchPre* const batchPre;
    precomp_x_t pre;

    // Replaces NewEvalPoint
    PreX(const BatchPre* const b, const fmpz_t x) : batchPre(b) {
        precomp_x_init(&pre, &batchPre->pre, x);
    }

    ~PreX() {
        precomp_x_clear(&pre);
    }

    void Eval(const fmpz_t* const yValues, fmpz_t out) {
        precomp_x_eval(&pre, yValues, out);
    }
};

struct CheckerPreComp {
    fmpz_t x;

    const BatchPre* degN;
    const BatchPre* deg2N;

    PreX *xN = nullptr;
    PreX *x2N = nullptr;

    CheckerPreComp(const size_t N)
    : degN(new BatchPre(roots, N))
    , deg2N(new BatchPre(roots2, 2 * N))
    {
        fmpz_init(x);
    }

    void setCheckerPrecomp(const fmpz_t val) {
        fmpz_set(x, val);

        clear_preX();
        xN = new PreX(degN, x);
        x2N = new PreX(deg2N, x);
    }

    void clear_preX() {
        if (xN != nullptr) {
            delete xN;
            delete x2N;
        }
        xN = nullptr;
        x2N = nullptr;
    }

    ~CheckerPreComp() {
        clear_preX();
        delete degN;
        delete deg2N;
        fmpz_clear(x);
    }
};

// TODO: Fix Both servers need to use same random seed. Using constant 1 instead now.
void randSum(fmpz_t out, const fmpz_t* const arr, const size_t len) {
    fmpz_t tmp;
    fmpz_init(tmp);

    for (unsigned int i = 0; i < len; i++) {
        // fmpz_randm(tmp, seed, Int_Modulus);
        fmpz_set_ui(tmp, 1);
        fmpz_mul(tmp, tmp, arr[i]);
        fmpz_mod(tmp, tmp, Int_Modulus);

        fmpz_add(out, out, tmp);
    }

    fmpz_mod(out, out, Int_Modulus);

    fmpz_clear(tmp);
}

struct Checker {
    const int server_num;     // id of this server
    const ClientPacket* const req;  // Client packet
    Circuit* const ckt;      // Validation circuit

    const size_t n;  // number of mult gates
    const size_t N;  // NextPowerofTwo(n)

    fmpz_t *pointsF;  // Points on f. f(i) = ith mul gate left input
    fmpz_t *pointsG;  // Points on g. g(i) = ith mul gate right input
    fmpz_t *pointsH;  // points on h. Want to check if h = f * g

    // For sigma = [r * (f(r) * g(r) - h(r))]
    fmpz_t evalF;  // [f(r)]
    fmpz_t evalG;  // [r * g(r)]
    fmpz_t evalH;  // [r * h(r)]

    Checker(Circuit* const c, const int idx, const ClientPacket* const req)
    : server_num(idx)
    , req(req)
    , ckt(c)
    , n(c->NumMulGates())
    , N(c->N())
    {
        new_fmpz_array(&pointsF, N);
        new_fmpz_array(&pointsG, N);
        new_fmpz_array(&pointsH, 2 * N);

        fmpz_init(evalF);
        fmpz_init(evalG);
        fmpz_init(evalH);

        ckt->ImportWires(req, server_num);
    }

    ~Checker() {
        clear_fmpz_array(pointsF, N);
        clear_fmpz_array(pointsG, N);
        clear_fmpz_array(pointsH, 2 * N);
        fmpz_clear(evalF);
        fmpz_clear(evalG);
        fmpz_clear(evalH);
    }

    void evalPoly(const CheckerPreComp* const pre) {
        // std::cout << "evalPoly" << std::endl;
        std::vector<Gate*> mulgates = ckt->mul_gates;
        // Get constant terms from packet
        fmpz_set(pointsF[0], req->f0_s);
        fmpz_set(pointsG[0], req->g0_s);
        fmpz_set(pointsH[0], req->h0_s);

        // For all multiplication triples a_i * b_i = c_i
        //    polynomial [f(x)] has [f(i)] = [a_i]
        //    polynomial [g(x)] has [g(i)] = [b_i]
        for (unsigned int i = 0; i < n; i++) {
            fmpz_set(pointsF[i+1], mulgates[i]->ParentL->WireValue);
            fmpz_set(pointsG[i+1], mulgates[i]->ParentR->WireValue);
            // Set even values of h to be output wires.
            fmpz_set(pointsH[2*(i + 1)], mulgates[i]->WireValue);
        }

        // Grab odd values of h from the packet.
        for (unsigned int j = 0; j < N; j++) {
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

    CorShare* CorShareFn(const CheckerPreComp* const pre) {
        evalPoly(pre);
        // std::cout << "CorShareFn" << std::endl;
        CorShare* out = new CorShare();

        fmpz_sub(out->shareD, evalF, req->triple_share->shareA);
        fmpz_mod(out->shareD, out->shareD, Int_Modulus);

        fmpz_sub(out->shareE, evalG, req->triple_share->shareB);
        fmpz_mod(out->shareE, out->shareE, Int_Modulus);

        return out;
    }

    void OutShare(fmpz_t out, const Cor* const corIn) const {
        fmpz_t mulCheck;
        fmpz_t term;

        fmpz_init(mulCheck);
        fmpz_init(term);

        if (server_num == 0) {
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
        const size_t num_zero_gates = ckt->result_zero.size();
        new_fmpz_array(&arr, num_zero_gates+1);

        fmpz_set(arr[0], mulCheck);

        for (unsigned int i = 0; i < num_zero_gates; i++) {
            fmpz_set(arr[i+1], ckt->result_zero[i]->WireValue);
        }

        randSum(out, arr, ckt->result_zero.size() + 1);

        clear_fmpz_array(arr, num_zero_gates+1);
        fmpz_clear(mulCheck);
        fmpz_clear(term);
    }
};

bool AddToZero(const fmpz_t x, const fmpz_t y) {
    fmpz_t sum; fmpz_init(sum);
    fmpz_add(sum, x, y);
    fmpz_mod(sum, sum, Int_Modulus);
    bool ans = fmpz_is_zero(sum);
    fmpz_clear(sum);
    return ans;
}

#endif
