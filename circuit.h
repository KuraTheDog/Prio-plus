#ifndef CIRCUIT_H
#define CIRCUIT_H

#include <iostream>
#include <vector>

#include "constants.h"
#include "share.h"

extern "C" {
    #include "poly/fft.h"
}

enum GateType {
    Gate_Input,
    Gate_Add,
    Gate_AddConst,
    Gate_Mul,
    Gate_MulConst
};

struct Gate {
    const GateType type;
    Gate* const ParentL;
    Gate* const ParentR;
    fmpz_t Constant;
    fmpz_t WireValue;

    Gate(GateType gatetype, Gate* pL = nullptr, Gate* pR = nullptr)
    : type(gatetype)
    , ParentL(pL)
    , ParentR(pR)
    {
        fmpz_init(Constant);
        fmpz_init(WireValue);
    }

    ~Gate() {
        fmpz_clear(Constant);
        fmpz_clear(WireValue);
    }
};

unsigned int NextPowerofTwo(const unsigned int n) {
    unsigned int ans = 1;
    while(n > ans)
        ans *= 2;
    return ans;
}

struct Circuit {
    std::vector<Gate*> gates;        // All gates
    // std::vector<Gate*> outputs;      // only used by unused
    std::vector<Gate*> result_zero;  // Gates that must be zero for eval to pass.
    // const size_t max_bits;           // Unused?

    std::vector<Gate*> mul_gates;
    std::vector<Gate*> mul_inp_gates;

    // Circuit(int n = 31) : max_bits(n) {}
    Circuit() {}

    ~Circuit() {
        // All gates, so don't need to also go over outputs and result_zero
        for (Gate* gate : gates)
            delete gate;
    }

    void addGate(Gate* gate) {
        gates.push_back(gate);
        if (gate->type == Gate_Mul) {
            mul_gates.push_back(gate);
            mul_inp_gates.push_back(gate);
        } else if (gate->type == Gate_Input) {
            mul_inp_gates.push_back(gate);
        }
    }

    void addZeroGate(Gate* zerogate) {
        result_zero.push_back(zerogate);
    }

    // Evals circuit on the input, returns if all result_zero gates are zero.
    bool Eval(const fmpz_t* const inps) {
        int inp_count = 0;
        for (unsigned int i = 0; i < gates.size(); i++) {

            switch (gates[i]->type)
            {
            case Gate_Input:
                fmpz_set(gates[i]->WireValue, inps[inp_count]);
                inp_count++;
                // std::cout << "  EVAL Input \tGate " << i << "  -  ";
                // fmpz_print(gates[i]->WireValue); std::cout << std::endl;
                break;
            case Gate_Add:
                fmpz_add(gates[i]->WireValue, gates[i]->ParentL->WireValue, gates[i]->ParentR->WireValue);
                fmpz_mod(gates[i]->WireValue, gates[i]->WireValue, Int_Modulus);
                // std::cout << "  EVAL Add \tGate " << i << "  -  ";
                // fmpz_print(gates[i]->ParentL->WireValue); std::cout << ", ";
                // fmpz_print(gates[i]->ParentR->WireValue); std::cout << " -> ";
                // fmpz_print(gates[i]->WireValue); std::cout << std::endl;
                break;
            case Gate_Mul:
                fmpz_mul(gates[i]->WireValue, gates[i]->ParentL->WireValue, gates[i]->ParentR->WireValue);
                fmpz_mod(gates[i]->WireValue, gates[i]->WireValue, Int_Modulus);
                // std::cout << "  EVAL Mul \tGate " << i << "  -  ";
                // fmpz_print(gates[i]->ParentL->WireValue); std::cout << ", ";
                // fmpz_print(gates[i]->ParentR->WireValue); std::cout << " -> ";
                // fmpz_print(gates[i]->WireValue); std::cout << std::endl;
                break;
            case Gate_AddConst:
                fmpz_add(gates[i]->WireValue, gates[i]->ParentL->WireValue, gates[i]->Constant);
                fmpz_mod(gates[i]->WireValue, gates[i]->WireValue, Int_Modulus);
                // std::cout << "  EVAL add con \tGate " << i << "  -  ";
                // fmpz_print(gates[i]->ParentL->WireValue); std::cout << ", const ";
                // fmpz_print(gates[i]->Constant); std::cout << " -> ";
                // fmpz_print(gates[i]->WireValue); std::cout << std::endl;
                break;
            case Gate_MulConst:
                fmpz_mul(gates[i]->WireValue, gates[i]->ParentL->WireValue, gates[i]->Constant);
                fmpz_mod(gates[i]->WireValue, gates[i]->WireValue, Int_Modulus);
                // std::cout << "  EVAL mul con \tGate " << i << "  -  ";
                // fmpz_print(gates[i]->ParentL->WireValue); std::cout << ", const ";
                // fmpz_print(gates[i]->Constant); std::cout << " -> ";
                // fmpz_print(gates[i]->WireValue); std::cout << std::endl;
                break;
            default:
                break;
            }
        }
        
        // all result_zero should be zero.
        for (Gate* zero_gate : this->result_zero) 
            if (not fmpz_is_zero(zero_gate->WireValue))
                return false;
        return true;
    }

    unsigned int NumMulGates() const {
        return mul_gates.size();
    }

    unsigned int N() const {
        return NextPowerofTwo(NumMulGates() + 1);
    }

    unsigned int NumMulInpGates() const {
        return mul_inp_gates.size();
    }

    void GetWireShares(fmpz_t** shares0, fmpz_t** shares1) const {
        unsigned int i = 0;

        for (Gate* gate : mul_inp_gates) {
            SplitShare(gate->WireValue, (*shares0)[i], (*shares1)[i]);
            i++;
        }
    }

    void ImportWires(const ClientPacket* const p, const int server_num) {
        unsigned int i = 0;

        for (Gate* gate : gates) {
            switch (gate->type) {
            case Gate_Input:
                fmpz_set(gate->WireValue, p->WireShares[i]);
                i++;
                break;
            case Gate_Add:
                fmpz_add(gate->WireValue, gate->ParentL->WireValue, gate->ParentR->WireValue);
                fmpz_mod(gate->WireValue, gate->WireValue, Int_Modulus);
                break;
            case Gate_Mul:
                fmpz_set(gate->WireValue, p->WireShares[i]);
                i++;
                break;
            case Gate_AddConst:
                if (server_num == 0)
                    fmpz_add(gate->WireValue, gate->ParentL->WireValue, gate->Constant);
                else
                    fmpz_set(gate->WireValue, gate->ParentL->WireValue);
                fmpz_mod(gate->WireValue, gate->WireValue, Int_Modulus);
                break;
            case Gate_MulConst:
                fmpz_mul(gate->WireValue, gate->ParentL->WireValue, gate->Constant);
                fmpz_mod(gate->WireValue, gate->WireValue, Int_Modulus);
                break;
            default:
                break;
            }
        }
    }
};

// Unused
/*
Circuit* AndCircuits(std::vector<Circuit*>& circuits) {
    Circuit* out = new Circuit();

    for (unsigned int i = 0; i < circuits.size(); i++) {
        out->gates.insert(out->gates.end(), circuits[i]->gates.begin(), circuits[i]->gates.end());
        out->outputs.insert(out->outputs.end(), circuits[i]->outputs.begin(), circuits[i]->outputs.end());
        out->result_zero.insert(out->result_zero.end(), circuits[i]->result_zero.begin(), circuits[i]->result_zero.end());
    }

    return out;
}
*/

Gate* MulByNegOne(Gate* const gate) {
    Gate* out = new Gate(Gate_MulConst, gate);

    fmpz_set_si(out->Constant, -1);
    fmpz_mod(out->Constant, out->Constant, Int_Modulus);

    return out;
}

/*
Various VALID(*) circuits.
*/

// Returns circuit that checks L*R == Prod
/*
Circuit* CheckMul(Gate* const L, Gate* const R, Gate* const Prod) {
    Circuit* out = new Circuit();

    Gate* mul = new Gate(Gate_Mul, L, R);

    Gate* inv = MulByNegOne(Prod);
    Gate* add = new Gate(Gate_Add, inv, mul);

    out->addGate(mul);
    out->addGate(inv);
    out->addGate(add);
    out->addZeroGate(add);

    return out;
}
*/

// Returns circuit for x^2 == y. For Varience and StdDev.
Circuit* CheckVar() {
    Gate* x = new Gate(Gate_Input);
    Gate* y = new Gate(Gate_Input);

    Circuit* out = new Circuit();

    out->addGate(x);
    out->addGate(y);

    Gate* mul = new Gate(Gate_Mul, x, x);

    Gate* inv = MulByNegOne(y);
    Gate* add = new Gate(Gate_Add, inv, mul);

    out->addGate(mul);
    out->addGate(inv);
    out->addGate(add);
    out->addZeroGate(add);

    return out;
}

/*
2: x, x^2, y, xy
3: x, x^2, x^3, x^4, y, xy, x^2 y
*/
Circuit* CheckLinReg(const size_t degree) {
    // std::cout << "Lin reg deg circuit: " << degree << std::endl;
    Circuit* out = new Circuit();
    
    const unsigned int num_inputs = 3 * (degree - 1) + 1;

    for (unsigned int i = 0; i < num_inputs; i++) {
        Gate* inp = new Gate(Gate_Input);
        out->addGate(inp);
    }

    Gate* x = out->gates[0];

    // Ensure x * x^i = x^{i+1}, or x * (x^i y) = x^{i+1} y
    for (unsigned int i = 0; i < 3 * (degree - 1); i++) {
        if (i == 2 * (degree - 1) - 1)  // skip x^d -> y
            continue;
        Gate* base = out->gates[i];
        Gate* next = out->gates[i+1];

        Gate* mul = new Gate(Gate_Mul, x, base);
        Gate* inv = MulByNegOne(next);
        Gate* add = new Gate(Gate_Add, mul, inv);

        out->addGate(mul);
        out->addGate(inv);
        out->addGate(add);
        out->addZeroGate(add);

        // std::cout << "Gate: inp[" << 0 << "] * inp[" << i << "] = inp[" << i+1 << "]" << std::endl;
    }

    return out;
}

#endif
