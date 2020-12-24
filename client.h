#ifndef CLIENT_H
#define CLIENT_H

#include "circuit.h"
#include "fmpz_utils.cpp"
#include "poly/fft.h"
#include "prio.h"
#include "share.h"

void share_polynomials(Circuit* circuit, ClientPacket& p0, ClientPacket& p1){
    auto mulgates = circuit->MulGates();

    int n = mulgates.size() + 1;

    int N = NextPowerofTwo(n);

    fmpz_t pointsF[N], pointsG[N];

    for(int i = 0; i < N; i++){
        fmpz_init(pointsF[i]);
        fmpz_init(pointsG[i]);
    }

    fmpz_t h0;
    fmpz_init(h0);

    fmpz_randm(pointsF[0],seed,Int_Modulus);
    fmpz_randm(pointsG[0],seed,Int_Modulus);
    fmpz_mul(h0,pointsF[0],pointsG[0]);
    fmpz_mod(h0,h0,Int_Modulus);

    for(int i = 1; i < n; i++){
        fmpz_set(pointsF[i],mulgates[i-1]->ParentL->WireValue);
        fmpz_set(pointsG[i],mulgates[i-1]->ParentR->WireValue);
    }

    for(int i = n; i < N; i++){
        fmpz_zero(pointsF[i]);
        fmpz_zero(pointsG[i]);
    }

    if(roots == nullptr){
        init_roots(N);
    }
    

    fmpz_t *polyF = fft_interpolate(Int_Modulus,N,invroots,pointsF,true);
    fmpz_t *polyG = fft_interpolate(Int_Modulus,N,invroots,pointsG,true);

    

    fmpz_t *paddedF, *paddedG;

    new_fmpz_array(&paddedF,2*N);
    new_fmpz_array(&paddedG,2*N);

    copy_fmpz_array(paddedF, polyF, N);
    copy_fmpz_array(paddedG, polyG, N);

    clear_fmpz_array(polyF, N);
    clear_fmpz_array(polyG, N);

    fmpz_t *evalsF = fft_interpolate(Int_Modulus,2*N,roots,paddedF,false);
    fmpz_t *evalsG = fft_interpolate(Int_Modulus,2*N,roots,paddedG,false);

    fmpz_t* h_points;
    new_fmpz_array(&h_points,N);
    
    int j = 0;
    std::cout << "hello" << std::endl;
    init_client_packet(p0, N);
    init_client_packet(p1, N);

    std::cout << "hello" << std::endl;

    for(int i = 1; i < 2*N-1; i+=2){
        fmpz_mul(h_points[j],evalsF[i],evalsG[i]);
        fmpz_mod(h_points[j],h_points[j], Int_Modulus);
        SplitShare(h_points[j],p0->h_points[j],p1->h_points[j]);
        j++;
    }
    

    SplitShare(pointsF[0],p0->f0_s,p1->f0_s);
    SplitShare(pointsG[0],p0->g0_s,p1->g0_s);
    SplitShare(h0,p0->h0_s,p1->h0_s);
    
    int num_wire_shares = 0;
    circuit->GetWireShares(&p0->WireShares,&p1->WireShares,num_wire_shares);
    
    clear_fmpz_array(pointsF,N);
    clear_fmpz_array(pointsG,N);
}
#endif