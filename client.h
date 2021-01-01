#ifndef CLIENT_H
#define CLIENT_H

#include "circuit.h"
#include "fmpz_utils.h"
#include "prio.h"
#include "share.h"

extern "C" {
    #include "poly/fft.h"
}

void share_polynomials(Circuit* circuit, ClientPacket& p0, ClientPacket& p1){
    auto mulgates = circuit->MulGates();

    int n = mulgates.size() + 1;

    int N = NextPowerofTwo(n);
    int NumMulInpGates = circuit->NumMulInpGates();

    // u_t, v_t = left and right wires of mul gates.
    // want f(t) = u_t, g(t) = v(t)
    fmpz_t* pointsF;
    fmpz_t* pointsG;
    new_fmpz_array(&pointsF, N);
    new_fmpz_array(&pointsG, N);

    fmpz_t h0;
    fmpz_init(h0);

    // Random f(0) = u_0, g(0) = v_0.
    fmpz_randm(pointsF[0],seed,Int_Modulus);
    fmpz_randm(pointsG[0],seed,Int_Modulus);
    // h(0) = f(0) * g(0)
    fmpz_mul(h0,pointsF[0],pointsG[0]);
    fmpz_mod(h0,h0,Int_Modulus);

    // u_i, v_i = left, right of i^th mult gate.
    for(int i = 1; i < n; i++){
        fmpz_set(pointsF[i],mulgates[i-1]->ParentL->WireValue);
        fmpz_set(pointsG[i],mulgates[i-1]->ParentR->WireValue);
    }

    // std::cout << " pointsF = [";
    // for (int i = 0; i < N; i++) {
    //     if (i > 0)
    //         std::cout << ", ";
    //     fmpz_print(pointsF[i]);
    // }
    // std::cout << "]" << std::endl;
    // std::cout << " pointsG = [";
    // for (int i = 0; i < N; i++) {
    //     if (i > 0)
    //         std::cout << ", ";
    //     fmpz_print(pointsG[i]);
    // }
    // std::cout << "]" << std::endl;

    // Initialize roots (nth roots of unity) and invroots (their inverse)
    if(roots == nullptr){
        init_roots(N);
    }

    // Build f, g that goes through pointsF, pointsG.
    // Interpolate through the Nth roots of unity
    fmpz_t *polyF = fft_interpolate(Int_Modulus,N,invroots,pointsF,true);
    fmpz_t *polyG = fft_interpolate(Int_Modulus,N,invroots,pointsG,true);

    // std::cout << " polyF = [";
    // for (int i = 0; i < N; i++) {
    //     if (i > 0) std::cout << ", ";
    //     fmpz_print(polyF[i]);
    // }
    // std::cout << "]" << std::endl;
    // std::cout << " polyG = [";
    // for (int i = 0; i < N; i++) {
    //     if (i > 0) std::cout << ", ";
    //     fmpz_print(polyG[i]);
    // }
    // std::cout << "]" << std::endl;

    // Pad to length 2N, to ensure it fits h.
    fmpz_t *paddedF, *paddedG;

    new_fmpz_array(&paddedF,2*N);
    new_fmpz_array(&paddedG,2*N);

    copy_fmpz_array(paddedF, polyF, N);
    copy_fmpz_array(paddedG, polyG, N);

    clear_fmpz_array(polyF, N);
    clear_fmpz_array(polyG, N);

    // Evaluate at all 2Nth roots of unity.
    fmpz_t *evalsF = fft_interpolate(Int_Modulus,2*N,roots2,paddedF,false);
    fmpz_t *evalsG = fft_interpolate(Int_Modulus,2*N,roots2,paddedG,false);

    // std::cout << " evalsF = [";
    // for (int i = 0; i < 2 * N; i++) {
    //     if (i > 0) std::cout << ", ";
    //     fmpz_print(evalsF[i]);
    // }
    // std::cout << "]" << std::endl;

    // std::cout << " evalsG = [";
    // for (int i = 0; i < 2 * N; i++) {
    //     if (i > 0) std::cout << ", ";
    //     fmpz_print(evalsG[i]);
    // }
    // std::cout << "]" << std::endl;
    
    init_client_packet(p0, N, NumMulInpGates);
    init_client_packet(p1, N, NumMulInpGates);

    // Send evaluations of f(r) * g(r) for all 2N-th roots of unity
    //     that aren't also N-th roots of unity
    // h_points[j] = evalF(2j + 1) * evalG(2j + 1), split into shares
    fmpz_t h_val;
    fmpz_init(h_val);
    for(int j = 0; j < N; j++){
        fmpz_mul(h_val, evalsF[2 * j + 1], evalsG[2 * j + 1]);
        fmpz_mod(h_val, h_val, Int_Modulus);
        SplitShare(h_val, p0->h_points[j], p1->h_points[j]);
    }
    
    // split f(0), g(0), h(0) into shares.
    SplitShare(pointsF[0],p0->f0_s,p1->f0_s);
    SplitShare(pointsG[0],p0->g0_s,p1->g0_s);
    SplitShare(h0,p0->h0_s,p1->h0_s);
    
    // Split outputs of input/mult gate shares.
    int num_wire_shares = 0;
    circuit->GetWireShares(&p0->WireShares,&p1->WireShares,num_wire_shares);
    
    auto triple = NewBeaverTriple();
    auto triple_shares = BeaverTripleShares(triple);

    p0->triple_share = &triple_shares[0];
    p1->triple_share = &triple_shares[1];

    delete triple;
    fmpz_clear(h0);
    fmpz_clear(h_val);
    clear_fmpz_array(pointsF,N);
    clear_fmpz_array(pointsG,N);
    clear_fmpz_array(paddedF, N);
    clear_fmpz_array(paddedG, N);
    clear_fmpz_array(evalsF, 2*N);
    clear_fmpz_array(evalsG, 2*N);
}

#endif
