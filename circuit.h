#ifndef CIRCUIT_H
#define CIRCUIT_H

#include <iostream>
#include <vector>

#include "constants.h"
#include "share.h"

extern "C" {
    #include "flint/fmpq_mat.h"  // for LinReg solving
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
    fmpz_t WireValue;  // set by eval/input, so gates not const

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

Gate* const MulByNegOne(Gate* const gate) {
    Gate* const out = new Gate(Gate_MulConst, gate);

    fmpz_set_si(out->Constant, -1);
    fmpz_mod(out->Constant, out->Constant, Int_Modulus);

    return out;
}

struct Circuit {
    std::vector<Gate*> gates;        // All gates
    // std::vector<Gate*> outputs;      // only used by unused
    std::vector<Gate*> result_zero;  // Gates that must be zero for eval to pass.
    // const size_t max_bits;           // Unused?

    std::vector<Gate*> mul_gates;

    // Circuit(int n = 31) : max_bits(n) {}
    Circuit() {}

    Circuit(const size_t num_inputs) : Circuit() {
        for (unsigned int i = 0; i < num_inputs; i++) {
            Gate* const inp = new Gate(Gate_Input);
            addGate(inp);
        }
    }

    ~Circuit() {
        // All gates, so don't need to also go over outputs and result_zero
        for (Gate* gate : gates)
            delete gate;
    }

    void addGate(Gate* const gate) {
        gates.push_back(gate);
        if (gate->type == Gate_Mul)
            mul_gates.push_back(gate);
    }

    void addZeroGate(Gate* const zerogate) {
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

    void GetMulShares(fmpz_t* const * const shares0,
                      fmpz_t* const * const shares1) const {
        unsigned int i = 0;

        for (Gate* gate : mul_gates) {
            SplitShare(gate->WireValue, (*shares0)[i], (*shares1)[i]);
            i++;
        }
    }

    void ImportWires(const ClientPacket* const p, const int server_num,
                     const fmpz_t* const InputShares) {
        unsigned int mul_idx = 0, inp_idx = 0;

        for (Gate* gate : gates) {
            switch (gate->type) {
            case Gate_Input:
                fmpz_set(gate->WireValue, InputShares[inp_idx]);
                inp_idx++;
                break;
            case Gate_Add:
                fmpz_add(gate->WireValue, gate->ParentL->WireValue, gate->ParentR->WireValue);
                fmpz_mod(gate->WireValue, gate->WireValue, Int_Modulus);
                break;
            case Gate_Mul:
                fmpz_set(gate->WireValue, p->MulShares[mul_idx]);
                mul_idx++;
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

    // Adds Gate[i] * Gate[j] = Gate[k] to circuit
    void AddCheckMulEqual(const size_t i, const size_t j, const size_t k) {
        Gate* const mul = new Gate(Gate_Mul, gates[i], gates[j]);
        Gate* const inv = MulByNegOne(gates[k]);
        Gate* const add = new Gate(Gate_Add, mul, inv);

        addGate(mul);
        addGate(inv);
        addGate(add);
        addZeroGate(add);

        // std::cout << "[" << i << "] * [" << j << "] = [" << k << "]\n";
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

/*
Various VALID(*) circuits.
*/

// Returns circuit for x^2 == y. For Varience and StdDev.
Circuit* const CheckVar() {
    Circuit* const out = new Circuit(2);

    out->AddCheckMulEqual(0, 0, 1);

    return out;
}

/*
2: x, x^2, y, xy
3: x0, x1, x0^2, x0x1, x1^2, y, x0 y, x1 y
4: x0, x1, x2, x0^2, x0x1, x0x2, x1^2, x1x2, x2^2, y, x0 y, x1 y, x2 y
*/
Circuit* const CheckLinReg(const size_t degree) {
    const size_t num_x = degree - 1;
    const size_t num_fields = 2 * num_x + 1 + num_x * (num_x + 1) / 2;

    Circuit* const out = new Circuit(num_fields);

    // xi * xj
    unsigned int idx = degree;
    for (unsigned int i = 0; i < num_x; i++) {
        for (unsigned int j = i; j < num_x; j++) {
            out->AddCheckMulEqual(i, j, idx);
            idx++;
        }
    }

    // y * xi
    for (unsigned int i = 0; i < num_x; i++)
        out->AddCheckMulEqual(num_x, i, idx + i);

    return out;
}

const double* const SolveLinReg(
        const size_t degree, const uint64_t* const x, const uint64_t* const y) {

    size_t idx = degree;

    // const size_t num_x = degree - 1;
    // std::cout << "n = " << x[0] << std::endl;
    // for (unsigned int i = 0; i < num_x; i++)
    //     std::cout << "x_" << i << " = " << x[i + 1] << std::endl;
    // std::cout << "y = " << y[0] << std::endl;

    // for (unsigned int i = 0; i < num_x; i++) {
    //     for (unsigned int j = i; j < num_x; j++) {
    //         std::cout << "x_" << i << " * " << "x_" << j << " = " << x[idx] << std::endl;
    //         idx++;
    //     }
    // }

    // for (unsigned int i = 0; i < num_x; i++)
    //     std::cout << "x_" << i << " * y = " << y[i + 1] << std::endl;

    fmpq_mat_t X; fmpq_mat_init(X, degree, degree);
    fmpq_mat_t Y; fmpq_mat_init(Y, degree, 1);

    for (unsigned int i = 0; i < degree; i++) {
        fmpq_set_ui(fmpq_mat_entry(X, i, 0), x[i], 1);
        if (i > 0)
            fmpq_set_ui(fmpq_mat_entry(X, 0, i), x[i], 1);

        fmpq_set_ui(fmpq_mat_entry(Y, i, 0), y[i], 1);
    }
    idx = degree;
    for (unsigned int i = 1; i < degree; i++) {
        for (unsigned int j = i; j < degree; j++) {
            fmpq_set_ui(fmpq_mat_entry(X, i, j), x[idx], 1);
            if (j != i)
                fmpq_set_ui(fmpq_mat_entry(X, j, i), x[idx], 1);
            idx++;
        }
    }

    // std::cout << "X: ";
    // fmpq_mat_print(X);
    // std::cout << "Y: ";
    // fmpq_mat_print(Y);

    int ret = fmpq_mat_inv(X, X);
    if (ret == 0) {
        std::cout << "WARNING: X is singular" << std::endl;
    }
    // std::cout << "X^-1: ";
    // fmpq_mat_print(X);

    fmpq_mat_mul(Y, X, Y);
    // std::cout << "X^-1 Y: ";
    // fmpq_mat_print(Y);


    double* ans = new double[degree];
    for (unsigned int i = 0; i < degree; i++)
        ans[i] = fmpq_get_d(fmpq_mat_entry(Y, i, 0));

    fmpq_mat_clear(X);
    fmpq_mat_clear(Y);

    return ans;
}

#endif
