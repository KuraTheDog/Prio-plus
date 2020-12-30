#ifndef SERVER_H
#define SERVER_H

#include "prio.h"
#include "share.h"
#include "circuit.h"
#include "fmpz_utils.h"

extern "C" {
    #include "poly/poly_batch.h"
    #include "poly/poly_once.h"
    #include "poly/fft.h"
}

struct BatchPoly {
    int nPoints;  // unused?
    fmpz_mod_poly_t fpoly;

    fmpz_t* Eval(fmpz_t* xPointsIn, int n) {
        std::cout << "Eval: n = " << n << std::endl << "xPointsIn = [";
        for (int i = 0; i < n; i++) {
            if (i > 0) std::cout << ", ";
            fmpz_print(xPointsIn[i]);
        }
        std::cout << "]" << std::endl;
        return poly_batch_evaluate(fpoly, n, xPointsIn);
    }

    void EvalOnce(fmpz_t x, fmpz_t out) {
        fmpz_init(out);
        poly_batch_evaluate_once(fpoly, x, out);
    }

    ~BatchPoly() {
        poly_batch_clear(fpoly);
    }
};

struct BatchPre {
    int nPoints;  // unused?
    precomp_t pre;

    // Newbatch
    BatchPre(fmpz_t *xPointsIn, int n) {
        nPoints = n;
        poly_batch_precomp_init(&pre, Int_Modulus, nPoints, xPointsIn);
    }

    // can we use this.npoints in place of n?
    BatchPoly* Interp(fmpz_t* yPointsIn, int n) {
        BatchPoly* bpoly = new BatchPoly();
        bpoly->nPoints = n;

        poly_batch_init(bpoly->fpoly, &pre);

        poly_batch_interpolate(bpoly->fpoly,&pre,yPointsIn);

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
    PreX(BatchPre* b, fmpz_t x) {
        batchPre = b;
        precomp_x_init(&pre, &batchPre->pre, x);
    }

    ~PreX() {
        precomp_x_clear(&pre);
    }

    void Eval(fmpz_t* yValues, fmpz_t out) {
        precomp_x_eval(&pre, yValues, out);
    }
};

struct CheckerPreComp {
    fmpz_t x;  // is also r?

    BatchPre *degN;
    BatchPre *deg2N;

    PreX *xN;
    PreX *x2N;

    CheckerPreComp(Circuit* ckt) {
        int n = ckt->NumMulGates() + 1;
        int N = NextPowerofTwo(n);

        degN = new BatchPre(roots, N);
        deg2N = new BatchPre(roots2, 2 * N);
    }

    void setCheckerPrecomp(fmpz_t val) {
        fmpz_init_set(x, val);

        xN = new PreX(degN, x);
        x2N = new PreX(deg2N, x);
    }
};

struct Checker {
    int serveridx;     // id of this server
    ClientPacket req;  // ??? Contains beaver triple?
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

    Checker(Circuit* c, int idx){
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

    void setReq(ClientPacket pkt){
        req = pkt;
        ckt->ImportWires(pkt, serveridx);
    }

    /* Step 2. 
       Given f0 and x can rebuild f. same with g. 
       Also rebuild wire shares.

       Can we replace pkt with req?
    */
    void evalPoly(CheckerPreComp *pre, ClientPacket pkt){
        std::vector<Gate*> mulgates = ckt->MulGates();
        // Get constant terms from packet
        fmpz_set(pointsF[0],pkt->f0_s);
        fmpz_set(pointsG[0],pkt->g0_s);
        fmpz_set(pointsH[0],pkt->h0_s);

        // For all multiplication triples a_i * b_i = c_i,
        //    polynomial [f(x)] has [f(i)] = [a_i]
        //    polynomial [g(x)] has [g(i)] = [b_i]
        for(int i = 1; i < n; i++){
            fmpz_set(pointsF[i],mulgates[i-1]->ParentL->WireValue);
            fmpz_set(pointsG[i],mulgates[i-1]->ParentR->WireValue);
            // Set even values of h to be output wires.
            fmpz_set(pointsH[2*i],mulgates[i-1]->WireValue);
        }

        // Can probably skip, since new_fmpz_array does this.
        for(int i = n; i < N; i++){
            fmpz_zero(pointsF[i]);
            fmpz_zero(pointsG[i]);
            fmpz_zero(pointsH[2*i]);
        }

        // Grab odd values of h from the packet.
        int j = 0;
        for(int i = 1; i < 2*N+1; i+=2){
            fmpz_set(pointsH[i], pkt->h_points[j]);
            j++;
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
};


#endif
