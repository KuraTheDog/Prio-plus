#ifndef SERVER_H
#define SERVER_H

#include "prio.h"
#include "share.h"
#include "circuit.h"

struct BatchPoly {
    int nPoints;
    fmpz_mod_poly_t fpoly;

    fmpz_t* Eval(fmpz_t* xPointsIn, int n){
        return poly_batch_evaluate(fpoly, n, xPointsIn);
    }

    fmpz_t EvalOnce(fmpz_t x){
        fmpz_t out;
        fmpz_init(out);
        poly_batch_evaluate_once(fpoly,x,out);
        return out;
    }
};

struct BatchPre {
    int nPoints;
    precomp_t pre;

    BatchPre(fmpz_t *xPointsIn, int n){
        nPoints = n;
        poly_batch_precomp_init(&pre, IntModulus, nPoints, xPointsIn);
    }

    BatchPoly* Interp(fmpz_t* yPointsIn, int n){
        BatchPoly* bpoly = (BatchPoly*) malloc(sizeof(BatchPoly));

        poly_batch_init(fpoly, &pre);

        poly_batch_interpolate(fpoly,&pre,yPointsIn);

        return bpoly;
    }

    
};

struct PreX {
    BatchPre *batchPre;
    precomp_x_t pre;
};

struct CheckerPreComp {
    fmpz_t x;

    BatchPre *degN;
    BatchPre *deg2N;

    PreX *xN;
    PreX *x2N;
};

struct Checker {
    int serveridx;
    ClientPacket req;
    Circuit *ckt;

    int n;
    int N;

    fmpz_t *pointsF;
    fmpz_t *pointsG;
    fmpz_t *pointsH;

    fmpz_t evalF;
    fmpz_t evalG;
    fmpz_t evalH;

    Checker(Circuit* c, int idx){
        ckt = c;
        serveridx = idx;
        n = ckt->gates.size()+1;
        N = NextPowerofTwo(n);

        pointsF = (fmpz_t *) malloc(N * sizeof(fmpz_t));
        pointsG = (fmpz_t *) malloc(N * sizeof(fmpz_t));
        pointsH = (fmpz_t *) malloc((2*N-1) * sizeof(fmpz_t));

        for(int i = 0; i < N; i++){
            fmpz_init(pointsF[i]);
            fmpz_init(pointsG[i]);
        }

        for(int i = 0; i < 2*N-1; i++){
            fmpz_init(pointsH[i]);
        }

        fmpz_init(evalF);
        fmpz_init(evalG);
        fmpz_init(evalH);
    }

    void setReq(ClientPacket pkt){
        req = pkt;
        ckt->ImportWires(pkt,serveridx);
    }

    void evalPoly(CheckerPreComp *pre, ClientPacket pkt){
        vector<Gate*> mulgates = ckt->MulGates();
        fmpz_set(pointsF[0],pkt->f0_s);
        fmpz_set(pointsG[0],pkt->g0_s);
        fmpz_set(pointsH[0],pkt->h0_s);

        for(int i = 1; i < n; i++){
            fmpz_set(pointsF[i],mulgates[i-1]->ParentL->WireValue);
            fmpz_set(pointsG[i],mulgates[i-1]->ParentR->WireValue);
            fmpz_set(pointsH[2*i],mulgates[i-1]->WireValue);
        }

        for(int i = n; i < N; i++){
            fmpz_zero(pointsF[i]);
            fmpz_zero(pointsG[i]);
            fmpz_zero(pointsH[2*i]);
        }

        int j = 0;
        for(int i = 1; i < 2*N-1; i+=2){
            fmpz_set(pointsH[i], pkt->h_points[j]);
            j++;
        }

        // set evals 
    }

};



#endif