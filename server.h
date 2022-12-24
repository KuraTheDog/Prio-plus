#ifndef SERVER_H
#define SERVER_H

#include "circuit.h"
#include "constants.h"
#include "interp.h"
#include "fmpz_utils.h"
#include "net_share.h"

// Return type of different ops
enum returnType {
    RET_INVALID,    // Too many inputs invalid
    RET_ANS,        // Success, Returning ans. For server returning an answer
    RET_NO_ANS,     // Success, no ans. For support server.
};

// Common (synced) randomness between all checkers on both servers
// Could also be static member, but then syncing is harder
extern flint_rand_t snips_seed;
flint_rand_t snips_seed;
void syncSnipSeeds(const int serverfd, const int server_num) {
    flint_randinit(snips_seed);
    fmpz_t tmp; fmpz_init(tmp);

    // fmpz_randm(tmp, snips_seed, Int_Modulus);
    // std::cout << "server " << server_num << " snips sync: next rand: ";
    // fmpz_print(tmp); std::cout << std::endl;

    if (server_num == 0)
        send_seed(serverfd, snips_seed);
    else
        recv_seed(serverfd, snips_seed);
}

struct Checker {
    const int server_num;     // id of this server
    const ClientPacket* const req;  // Client packet
    Circuit* const ckt;      // Validation circuit

    const size_t n;  // number of mult gates
    const size_t N;  // NextPowerOfTwo(n)

    fmpz_t* pointsF;  // Points on f. f(i) = ith mul gate left input
    fmpz_t* pointsG;  // Points on g. g(i) = ith mul gate right input
    fmpz_t* pointsH;  // points on h. Want to check if h = f * g

    // For sigma = [r * (f(r) * g(r) - h(r))]
    fmpz_t evalF;  // [f(r)]
    fmpz_t evalG;  // [r * g(r)]
    fmpz_t evalH;  // [r * h(r)]

    // if both servers on same runtime/execution (e.g. test files), use fixed randomness instead
    const bool same_runtime = false;

    Checker(Circuit* const c, const int idx, const ClientPacket* const req,
            const MultCheckPreComp* const pre, const fmpz_t* const InputShares,
            const bool same_runtime = false)
    : server_num(idx)
    , req(req)
    , ckt(c)
    , n(c->NumMulGates())
    , N(NextPowerOfTwo(n))
    , same_runtime(same_runtime)
    {
        if (same_runtime) {
            std::cout << "DEBUG: using fixed checker randomness since same runtime" << std::endl;
            flint_randinit(snips_seed);
        }

        new_fmpz_array(&pointsF, N + 1);
        new_fmpz_array(&pointsG, N + 1);
        new_fmpz_array(&pointsH, 2 * (N + 1));

        fmpz_init(evalF);
        fmpz_init(evalG);
        fmpz_init(evalH);

        ckt->ImportWires(req, server_num, InputShares);
        evalPoly(pre);
    }

    ~Checker() {
        flint_randclear(seed);
        clear_fmpz_array(pointsF, N + 1);
        clear_fmpz_array(pointsG, N + 1);
        clear_fmpz_array(pointsH, 2 * (N + 1));
        fmpz_clear(evalF);
        fmpz_clear(evalG);
        fmpz_clear(evalH);
    }

    void evalPoly(const MultCheckPreComp* const pre) {
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

        // set evals, extra times x to poly
        fmpz_t x; fmpz_init(x); pre->getEvalPoint(x);
        pre->Eval(pointsF, evalF);
        pre->Eval(pointsG, evalG);
        fmpz_mul(evalG, evalG, x);
        fmpz_mod(evalG, evalG, Int_Modulus);
        pre->Eval2(pointsH, evalH);
        fmpz_mul(evalH, evalH, x);
        fmpz_mod(evalH, evalH, Int_Modulus);
        fmpz_clear(x);
    }

    CorShare* CorShareFn() {
        // std::cout << "CorShareFn" << std::endl;
        CorShare* const out = new CorShare();

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

    void randSum(fmpz_t out, const fmpz_t* const arr, const size_t len) const {
        fmpz_t tmp; fmpz_init(tmp);

        for (unsigned int i = 0; i < len; i++) {
            fmpz_randm(tmp, snips_seed, Int_Modulus);
            if (same_runtime)
                fmpz_set_ui(tmp, 1);
            fmpz_mul(tmp, tmp, arr[i]);
            fmpz_mod(tmp, tmp, Int_Modulus);

            fmpz_add(out, out, tmp);
        }

        fmpz_mod(out, out, Int_Modulus);

        fmpz_clear(tmp);
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
