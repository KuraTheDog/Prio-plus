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
fmpz_t *roots = nullptr, *invroots = nullptr, *roots2 = nullptr;

// TODO: move these next 3 things to a prio.cpp for initializing the globals.

void init_constants() {
    fmpz_init(Int_Modulus);
    fmpz_set_str(Int_Modulus,Int_Modulus_str.c_str(),16);
    fmpz_init(Int_Gen);
    fmpz_set_str(Int_Gen,Int_Gen_str.c_str(),16);

    std::cout << "Init constants: " << std::endl;
    std::cout << "  Int_Modulus = "; fmpz_print(Int_Modulus); std::cout << std::endl;
    std::cout << "  Int_Gen = "; fmpz_print(Int_Gen); std::cout << std::endl;
    std::cout << "  twoOrder = " << twoOrder << std::endl;
    flint_randinit(seed);
}

void clear_constants() {
    flint_randclear(seed);
}

void init_roots(const int N) {
    // std::cout << "init_roots: " << N << std::endl;

    new_fmpz_array(&roots, N);
    new_fmpz_array(&invroots, N);
    new_fmpz_array(&roots2, 2 * N);

    int step_size = (1 << twoOrder)/N;  // 2^(twoOrder - log_2 N)
    fmpz_t g_, ginv_, ghalf_;
    fmpz_init(g_);      // Generator (Int_Gen)
    fmpz_init(ginv_);   // Inverse of g_
    fmpz_init(ghalf_);  // g_^(step/2), for 2N roots.

    fmpz_invmod(ginv_,Int_Gen,Int_Modulus);

    /*
    N = 2^k, so stepsize = 2^(Ord - k).
    g_ = gen^stepsize.
    So g_^N = gen^(2^ord) = 1, by fermat little. (Actually 2^ord - 1 enough?)
    */
    fmpz_powm_ui(g_,Int_Gen,step_size, Int_Modulus);
    fmpz_powm_ui(ginv_,ginv_,step_size, Int_Modulus);
    fmpz_powm_ui(ghalf_, Int_Gen, step_size / 2, Int_Modulus);
    fmpz_set_ui(roots[0],1);
    fmpz_set_ui(invroots[0],1);
    fmpz_set_ui(roots2[0],1);

    for(int i = 1; i < N; i++){
        fmpz_mul(roots[i],roots[i-1],g_);
        fmpz_mul(invroots[i],invroots[i-1],ginv_);

        fmpz_mod(roots[i],roots[i],Int_Modulus);
        fmpz_mod(invroots[i],invroots[i],Int_Modulus);
    }

    for(int i = 1; i < 2 * N; i++){
        fmpz_mul(roots2[i], roots2[i-1], ghalf_);
        fmpz_mod(roots2[i],roots2[i],Int_Modulus);
    }

    std::cout << " roots = {";
    for (int i = 0; i < N; i++) {
        if (i > 0) std::cout << ", ";
        fmpz_print(roots[i]);
    }
    std::cout << "}" << std::endl;

    std::cout << " invroots = {";
    for (int i = 0; i < N; i++) {
        if (i > 0) std::cout << ", ";
        fmpz_print(invroots[i]);
    }
    std::cout << "}" << std::endl;

    std::cout << " roots2 = {";
    for (int i = 0; i < 2 * N; i++) {
        if (i > 0) std::cout << ", ";
        fmpz_print(roots2[i]);
    }
    std::cout << "}" << std::endl;

    fmpz_clear(g_);
    fmpz_clear(ginv_);
    fmpz_clear(ghalf_);
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

    ~Gate() {
        fmpz_clear(Constant);
        fmpz_clear(WireValue);
    }
};

struct Circuit {
    std::vector<Gate*> gates;        // All gates
    std::vector<Gate*> outputs;      // ?
    std::vector<Gate*> result_zero;  // Gates that must be zero for eval to pass.
    int max_bits;                    // 

    Circuit(int n = 31){
        max_bits = n;
    }

    void addGate(Gate* gate){
        gates.push_back(gate);
    }

    void addZeroGate(Gate* zerogate){
        result_zero.push_back(zerogate);
    }

    // Evals circuit on the input, returns if all result_zero gates are zero.
    bool Eval(fmpz_t *inps){
        int inp_count = 0;
        // std::cout << "EVAL" << std::endl;
        for(int i = 0; i < gates.size(); i++){
            
            switch (gates[i]->type)
            {
            case Gate_Input:
                fmpz_set(gates[i]->WireValue,inps[inp_count]);
                inp_count++;
                // std::cout << "  EVAL Input \tGate " << i << "  -  ";
                // fmpz_print(gates[i]->WireValue); std::cout << std::endl;
                break;
            case Gate_Add:
                fmpz_add(gates[i]->WireValue,gates[i]->ParentL->WireValue,gates[i]->ParentR->WireValue);
                fmpz_mod(gates[i]->WireValue,gates[i]->WireValue,Int_Modulus);
                // std::cout << "  EVAL Add \tGate " << i << "  -  ";
                // fmpz_print(gates[i]->ParentL->WireValue); std::cout << ", ";
                // fmpz_print(gates[i]->ParentR->WireValue); std::cout << " -> ";
                // fmpz_print(gates[i]->WireValue); std::cout << std::endl;
                break;
            case Gate_Mul:
                fmpz_mul(gates[i]->WireValue,gates[i]->ParentL->WireValue,gates[i]->ParentR->WireValue);
                fmpz_mod(gates[i]->WireValue,gates[i]->WireValue,Int_Modulus);
                // std::cout << "  EVAL Mul \tGate " << i << "  -  ";
                // fmpz_print(gates[i]->ParentL->WireValue); std::cout << ", ";
                // fmpz_print(gates[i]->ParentR->WireValue); std::cout << " -> ";
                // fmpz_print(gates[i]->WireValue); std::cout << std::endl;
                break;
            case Gate_AddConst:
                fmpz_add(gates[i]->WireValue,gates[i]->ParentL->WireValue,gates[i]->Constant);
                fmpz_mod(gates[i]->WireValue,gates[i]->WireValue,Int_Modulus);
                // std::cout << "  EVAL add con \tGate " << i << "  -  ";
                // fmpz_print(gates[i]->ParentL->WireValue); std::cout << ", const ";
                // fmpz_print(gates[i]->Constant); std::cout << " -> ";
                // fmpz_print(gates[i]->WireValue); std::cout << std::endl;
                break;
            case Gate_MulConst:
                fmpz_mul(gates[i]->WireValue,gates[i]->ParentL->WireValue,gates[i]->Constant);
                fmpz_mod(gates[i]->WireValue,gates[i]->WireValue,Int_Modulus);
                // std::cout << "  EVAL mul con \tGate " << i << "  -  ";
                // fmpz_print(gates[i]->ParentL->WireValue); std::cout << ", const ";
                // fmpz_print(gates[i]->Constant); std::cout << " -> ";
                // fmpz_print(gates[i]->WireValue); std::cout << std::endl;
                break;
            default:
                break;
            }
        }
        bool res = true;

        for(auto zero_gate : this->result_zero)
            res = res and fmpz_is_zero(zero_gate->WireValue);

        return res;
    }

    std::vector<Gate*> MulGates() const {
        std::vector<Gate*> res;

        for(auto gate : this->gates){
            if(gate->type == Gate_Mul)
                res.push_back(gate);
        }

        return res;
    }

    int NumMulGates() const {
        int total = 0;
        for (auto gate: this->gates) {
            if (gate->type == Gate_Mul)
                total += 1;
        }
        return total;
    }

    int NumMulInpGates() const {
        int total = 0;
        for (auto gate: this->gates) {
            if (gate->type == Gate_Mul or gate->type == Gate_Input)
                total += 1;
        }
        return total;
    }

    void GetWireShares(fmpz_t** shares0, fmpz_t** shares1, int& n) const {
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

        // deal with input gates modification as they should be XOR shared for length checks
        for(auto gate : gates) {
            if(gate->type == Gate_Mul or gate->type == Gate_Input){
                SplitShare(gate->WireValue,forServer0[i],forServer1[i]);
                i++;
            }
        }

        *shares0 = forServer0; *shares1 = forServer1;
    }

    void ImportWires(const ClientPacket p, const int server_num){
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

int NextPowerofTwo(const int n){
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

    fmpz_set_si(out->Constant,-1);
    fmpz_mod(out->Constant,out->Constant,Int_Modulus);

    out->ParentL = gate;

    return out;
}

/*
Various VALID(*) circuits.
*/

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


// Returns circuit for x^2 == y. For Varience and StdDev.
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
