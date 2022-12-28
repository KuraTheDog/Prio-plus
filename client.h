#ifndef CLIENT_H
#define CLIENT_H

#include "circuit.h"
#include "constants.h"
#include "interp.h"
#include "fmpz_utils.h"
#include "share.h"

// Expects p0, p1 to already be initialized
void share_polynomials(const Circuit* const circuit, ClientPacket* const p0, ClientPacket* const p1) {
    auto mulgates = circuit->mul_gates;

    const unsigned int n = circuit->NumMulGates();
    const unsigned int N = NextPowerOfTwo(n);

    // u_t, v_t = left and right wires of mul gates.
    // want f(t) = u_t, g(t) = v(t)
    fmpz_t* pointsF; new_fmpz_array(&pointsF, N);
    fmpz_t* pointsG; new_fmpz_array(&pointsG, N);

    fmpz_t h0; fmpz_init(h0);

    // Random f(0) = u_0, g(0) = v_0.
    fmpz_randm(pointsF[0], seed, Int_Modulus);
    fmpz_randm(pointsG[0], seed, Int_Modulus);
    // h(0) = f(0) * g(0)
    fmpz_mul(h0, pointsF[0], pointsG[0]);
    fmpz_mod(h0, h0, Int_Modulus);

    // u_i, v_i = left, right of i^th mult gate.
    for (unsigned int i = 0; i < n; i++) {
        fmpz_set(pointsF[i + 1], mulgates[i]->ParentL->WireValue);
        fmpz_set(pointsG[i + 1], mulgates[i]->ParentR->WireValue);
    }

    // Build f, g that goes through pointsF, pointsG.
    // Interpolate through the Nth roots of unity
    fmpz_t* const polyF = interpolate_inv(N, pointsF);
    fmpz_t* const polyG = interpolate_inv(N, pointsG);

    // Pad to length 2N, to ensure it fits h.
    fmpz_t* paddedF; new_fmpz_array(&paddedF, 2*N);
    fmpz_t* paddedG; new_fmpz_array(&paddedG, 2*N);

    copy_fmpz_array(paddedF, polyF, N);
    copy_fmpz_array(paddedG, polyG, N);

    clear_fmpz_array(polyF, N);
    clear_fmpz_array(polyG, N);

    // Evaluate at all 2Nth roots of unity.
    fmpz_t* const evalsF = interpolate_2N(N, paddedF);
    fmpz_t* const evalsG = interpolate_2N(N, paddedG);

    // Send evaluations of f(r) * g(r) for all 2N-th roots of unity
    //     that aren't also N-th roots of unity
    // h_points[j] = evalF(2j + 1) * evalG(2j + 1), split into shares
    fmpz_t h_val;
    fmpz_init(h_val);
    for (unsigned int j = 0; j < N; j++) {
        fmpz_mul(h_val, evalsF[2 * j + 1], evalsG[2 * j + 1]);
        fmpz_mod(h_val, h_val, Int_Modulus);
        SplitShare(h_val, p0->h_points[j], p1->h_points[j]);
    }

    // split f(0), g(0), h(0) into shares.
    SplitShare(pointsF[0], p0->f0_s, p1->f0_s);
    SplitShare(pointsG[0], p0->g0_s, p1->g0_s);
    SplitShare(h0, p0->h0_s, p1->h0_s);

    // Split outputs of input/mult gate shares.
    circuit->GetMulShares(&p0->MulShares, &p1->MulShares);

    NewBeaverTriples(p0->triple, p1->triple);

    fmpz_clear(h0);
    fmpz_clear(h_val);
    clear_fmpz_array(pointsF, N);
    clear_fmpz_array(pointsG, N);
    clear_fmpz_array(paddedF, N);
    clear_fmpz_array(paddedG, N);
    clear_fmpz_array(evalsF, 2*N);
    clear_fmpz_array(evalsG, 2*N);
}

#endif
