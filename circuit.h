#ifndef CIRCUIT_H
#define CIRCUIT_H

#include "fmpz_utils.h"
#include "prio.h"
#include "share.h"

extern "C" {
    #include "poly/fft.h"
}

/*
Globals (extern) from prio.h,
Modulus and seed are also used (declared) in share.cpp
If they're here too, this causes redefinition issues.

Ideally, this would all be in some circuit.cpp, rather than header file, which would fix this issue from happening.
*/
// fmpz_t Int_Modulus;
fmpz_t Int_Gen;
// flint_rand_t seed;
fmpz_t *roots = nullptr, *invroots = nullptr;

// TODO: move these next 3 things to a prio.cpp for initializing the globals.

void init_constants() {
    fmpz_init(Int_Modulus);
    fmpz_set_str(Int_Modulus,Int_Modulus_str.c_str(),16);
    fmpz_init(Int_Gen);
    fmpz_set_str(Int_Gen,Int_Gen_str.c_str(),16);

    flint_randinit(seed);
}

void clear_constants() {
    flint_randclear(seed);
}

void init_roots(int N) {

    roots = (fmpz_t *) malloc(N* sizeof(fmpz_t));
    invroots = (fmpz_t *) malloc(N* sizeof(fmpz_t));

    init_fmpz_array(roots, N);
    init_fmpz_array(invroots, N);


    int step_size = (1 << twoOrder)/N;
    fmpz_t g_, ginv_;
    fmpz_init(g_);
    fmpz_init(ginv_);

    fmpz_invmod(ginv_,Int_Gen,Int_Modulus);

    fmpz_pow_ui(g_,Int_Gen,step_size);
    fmpz_pow_ui(ginv_,ginv_,step_size);
    fmpz_set_ui(roots[0],1);
    fmpz_set_ui(invroots[0],1);

    for(int i = 1; i < N; i++){
        fmpz_mul(roots[i],roots[i-1],g_);
        fmpz_mul(invroots[i],roots[i-1],ginv_);

        fmpz_mod(roots[i],roots[i],Int_Modulus);
        fmpz_mod(invroots[i],invroots[i],Int_Modulus);
    }

}

enum GateType {
    Gate_Input,
    Gate_Add,
    Gate_AddConst,
    Gate_Mul,
    Gate_MulConst
};

struct Gate {
    GateType type;
    Gate* ParentL;
    Gate* ParentR;
    fmpz_t Constant;
    fmpz_t WireValue;

    Gate(){
        ParentL = ParentR = nullptr;
        fmpz_init(Constant);
        fmpz_init(WireValue);
    }

    Gate(GateType gatetype) {
        type = gatetype;
        ParentL = ParentR = nullptr;
        fmpz_init(Constant);
        fmpz_init(WireValue);
    }

};

struct Circuit {
    std::vector<Gate*> gates;
    std::vector<Gate*> outputs;
    std::vector<Gate*> result_zero;

    Circuit(){
    }

    void addGate(Gate* gate){
        gates.push_back(gate);
    }

    void addZeroGate(Gate* zerogate){
        result_zero.push_back(zerogate);
    }

    bool Eval(fmpz_t *inps){
        int inp_count = 0;
        std::cout << "EVAL" << std::endl;
        for(int i = 0; i < gates.size(); i++){
            
            switch (gates[i]->type)
            {
            case Gate_Input:
                fmpz_set(gates[i]->WireValue,inps[inp_count]);
                inp_count++;
                break;
            case Gate_Add:
                fmpz_add(gates[i]->WireValue,gates[i]->ParentL->WireValue,gates[i]->ParentR->WireValue);
                fmpz_mod(gates[i]->WireValue,gates[i]->WireValue,Int_Modulus);
                break;
            case Gate_Mul:
                fmpz_mul(gates[i]->WireValue,gates[i]->ParentL->WireValue,gates[i]->ParentR->WireValue);
                fmpz_mod(gates[i]->WireValue,gates[i]->WireValue,Int_Modulus);
                break;
            case Gate_AddConst:
                fmpz_add(gates[i]->WireValue,gates[i]->ParentL->WireValue,gates[i]->Constant);
                fmpz_mod(gates[i]->WireValue,gates[i]->WireValue,Int_Modulus);
                break;
            case Gate_MulConst:
                fmpz_mul(gates[i]->WireValue,gates[i]->ParentL->WireValue,gates[i]->Constant);
                fmpz_mod(gates[i]->WireValue,gates[i]->WireValue,Int_Modulus);
                break;
            default:
                break;
            }
            std::cout << "EVAL Gate " << i << "  -  ";
            fmpz_print(gates[i]->WireValue);
            std::cout << std::endl;
        }
        bool res = true;

        for(auto zero_gate : this->result_zero)
            res = res and fmpz_is_zero(zero_gate->WireValue);

        return res;
    }

    std::vector<Gate*> MulGates(){
        std::vector<Gate*> res;

        for(auto gate : this->gates){
            if(gate->type == Gate_Mul)
                res.push_back(gate);
        }

        return res;
    }

    void GetWireShares(fmpz_t** shares0, fmpz_t** shares1, int& n){
        n = 0;
        for(auto gate : gates){
            if(gate->type == Gate_Input or gate->type == Gate_Mul){
                n++;
            }
        }

        fmpz_t* forServer0 = (fmpz_t *) malloc(n*sizeof(fmpz_t));
        fmpz_t* forServer1 = (fmpz_t *) malloc(n*sizeof(fmpz_t));

        for(int i = 0; i < n; i++){
            fmpz_init(forServer0[i]);
            fmpz_init(forServer1[i]);
        }

        int i = 0;

        for(auto gate : gates) {
            if(gate->type == Gate_Input or gate->type == Gate_Mul){
                SplitShare(gate->WireValue,forServer0[i],forServer1[i]);
                i++;
            }
        }

        *shares0 = forServer0; *shares1 = forServer1;
    }

    void ImportWires(ClientPacket p, int server_num){
        int n = 0;
        for(auto gate : gates){
            if(gate->type == Gate_Input or gate->type == Gate_Mul){
                n++;
            }
        }
        int i = 0;


        for(auto gate : gates){
            switch (gate->type)
            {
            case Gate_Input:
                fmpz_set(gate->WireValue,p->WireShares[i]);
                i++;
                break;
            case Gate_Add:
                fmpz_add(gate->WireValue,gate->ParentL->WireValue,gate->ParentR->WireValue);
                fmpz_mod(gate->WireValue,gate->WireValue,Int_Modulus);
                break;
            case Gate_Mul:
                fmpz_set(gate->WireValue,p->WireShares[i]);
                i++;
                break;
            case Gate_AddConst:
                if(server_num == 0)
                    fmpz_add(gate->WireValue,gate->ParentL->WireValue,gate->Constant);
                else
                    fmpz_set(gate->WireValue,gate->ParentL->WireValue);
                fmpz_mod(gate->WireValue,gate->WireValue,Int_Modulus);
                break;
            case Gate_MulConst:
                fmpz_mul(gate->WireValue,gate->ParentL->WireValue,gate->Constant);
                fmpz_mod(gate->WireValue,gate->WireValue,Int_Modulus);
                break;
            default:
                break;
            }
        }
    }
};

int NextPowerofTwo(int n){
    int ans = 1;
    while(n > ans){
        ans *= 2;
    }

    return ans;
}

Circuit* AndCircuits(std::vector<Circuit*>& circuits) {
    Circuit* out = new Circuit();

    for(int i = 0; i < circuits.size(); i++){
        out->gates.insert(out->gates.end(),circuits[i]->gates.begin(),circuits[i]->gates.end());
        out->outputs.insert(out->outputs.end(),circuits[i]->outputs.begin(),circuits[i]->outputs.end());
        out->result_zero.insert(out->result_zero.end(),circuits[i]->result_zero.begin(),circuits[i]->result_zero.end());
    }

    return out;
}

Gate* MulByNegOne(Gate* gate) {
    Gate* out = new Gate(Gate_MulConst);

    std::cout << "Negone" << std::endl;

    fmpz_set_si(out->Constant,-1);
    fmpz_mod(out->Constant,out->Constant,Int_Modulus);

    out->ParentL = gate;

    return out;
}

// Returns circuit that checks L*R == Prod
Circuit* CheckMul(Gate* L, Gate* R, Gate* Prod) {
    Circuit* out = new Circuit();

    Gate* mul = new Gate(Gate_Mul);
    mul->ParentL = L;
    mul->ParentR = R;

    Gate* inv = MulByNegOne(Prod);
    Gate* add = new Gate(Gate_Add);
    add->ParentL = inv;
    add->ParentR = mul;

    out->addGate(mul);
    out->addGate(inv);
    out->addGate(add);
    out->addZeroGate(add);

    return out;
}

Circuit* CheckVar(){
    Gate* x = new Gate(Gate_Input);
    Gate* y = new Gate(Gate_Input);

    Circuit* out = new Circuit();

    out->addGate(x);
    out->addGate(y);

    Gate* mul = new Gate(Gate_Mul);
    mul->ParentL = x;
    mul->ParentR = x;

    Gate* inv = MulByNegOne(y);
    Gate* add = new Gate(Gate_Add);
    add->ParentL = inv;
    add->ParentR = mul;

    out->addGate(mul);
    out->addGate(inv);
    out->addGate(add);
    out->addZeroGate(add);

    return out;
}



// Circuit* CheckLinReg() {

// }

#endif
