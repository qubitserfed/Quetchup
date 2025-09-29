#include <algorithm>
#include <complex>
#include <fstream>
#include <iostream>
#include <functional>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include <cmath>

#include "linear_algebra.hpp"
#include "clifford_state.hpp"
#include "clifford_map.hpp"

using f64 = double;
using Complex = std::complex<f64>;

const f64 PI = M_PI;
const Complex I = Complex(0, 1);

struct Gate {
    std::string name;
    int i, j;
};


std::vector<std::string> gate_types = {
    "H",
    "S",
    "CZ",
    "SWAP",
    "POP",
    "PUSH"
};

template<typename T>
T sample(std::vector<T> vec) {
    return vec[rand() % vec.size()];
}

std::vector<Gate> rand_circ(int qubits, int depth) {
    std::vector<Gate> res;
    int live_qubits = qubits;
    for (int i = 0; i < depth; ++i) {
        int a, b;
        std::string gate;
        if (live_qubits == 1) {
            std::vector<std::string> gate_types = { "H", "S", "X", "Z" };
            if (qubits > 1)
                gate_types.push_back("PUSH");
            a = 0;
            b = 0;
            gate = sample(gate_types);
            if (gate == "PUSH") {
                live_qubits++;
            }
        }
        else {
            std::vector<std::string> gate_types = { "CZ", "SWAP", "POP", "H", "S", "CZ" };
            if (live_qubits < qubits)
                gate_types.push_back("PUSH");
            gate = sample(gate_types);

            if (gate == "POP") {
                live_qubits--;
            }
            else if (gate == "PUSH") {
                live_qubits++;
            }
            b = rand() % (live_qubits - 1) + 1;
            a = rand() % b;
        }
        res.push_back({gate, a, b});
    }

    return res;
}

int clifford_group_size(int N) {
    std::vector<CliffordMap> group;

    std::function<CliffordMap(int)> h_i = [&](int idx) {
        CliffordMap res = id_map(N);
        apply_h(res.state, idx + N);
        return res;
    };

    std::function<CliffordMap(int)> s_i = [&](int idx) {
        CliffordMap res = id_map(N);
        apply_s(res.state, idx + N);
        return res;
    };
    
    std::function<CliffordMap(int, int)> cz_ij = [&](int idx0, int idx1) {
        CliffordMap res = id_map(N);
        apply_cz(res.state, idx0 + N, idx1 + N);
        return res;
    };
    
    std::set<std::string> group_strings;
    
    group.push_back(id_map(N));

    int layer_l = 0;
    int layer_r = 1;
    while (layer_l < layer_r) {
        // apply S gates
        std::cout << layer_r << std::endl;
        for (int i = 0; i < N; ++i) {
            CliffordMap si = s_i(i);
            for (int j = layer_l; j < layer_r; ++j) {
                CliffordMap gj = group[j];
                CliffordMap gj_si = compose(gj, si);
                if (group_strings.find(to_string(gj_si.state)) == group_strings.end()) {
                    group_strings.insert(to_string(gj_si.state));
                    group.push_back(gj_si);
                }
            }
        }

        // apply H gates
        for (int i = 0; i < N; ++i) {
            CliffordMap hi = h_i(i);
            for (int j = layer_l; j < layer_r; ++j) {
                CliffordMap gj = group[j];
                CliffordMap gj_hi = compose(gj, hi);
                if (group_strings.find(to_string(gj_hi.state)) == group_strings.end()) {
                    group_strings.insert(to_string(gj_hi.state));
                    group.push_back(gj_hi);
                }
            }
        }

        // apply CZ gates
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) if (i != j) {
                CliffordMap czij = cz_ij(i, j);
                for (int k = layer_l; k < layer_r; ++k) {
                    CliffordMap gk = group[k];
                    CliffordMap gk_czij= compose(gk, czij);
                    if (group_strings.find(to_string(gk_czij.state)) == group_strings.end()) {
                        group_strings.insert(to_string(gk_czij.state));
                        group.push_back(gk_czij);
                    }
                }
            }
        }

        layer_l = layer_r;
        layer_r = group.size();
        
    };


    int p[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    std::cout << '\n';
    for(int i = 0; i < group.size(); ++i) {
        p[group[i].state.phase]++;
        continue;
        auto g = group[i];
        std::cout << i << ": \n";
        print(g.state);
        print_superposition(g.state);
        std::cout << "--------------------------------" << std::endl;
    }
    for (int i = 0; i < 8; ++i)
        std::cout << p[i] << " ";
    std::cout << std::endl;

    return group.size();
}

struct BruteState {
    int n; // number of qubits
    std::vector<Complex> vals;
};

Complex to_root8(int k) {
    return Complex(cos(2.0 * PI * k / 8), sin(2.0 * PI * k / 8));
}

int from_root8(Complex cpx) {
    int res = -1;
    f64 mindist = std::numeric_limits<f64>::infinity();

    for (int k = 0; k < 8; ++k) {
        f64 dist = abs(cpx - to_root8(k));
        if (dist < mindist) {
            mindist = dist;
            res = k;
        }
    }
    return res;
}

void clean(BruteState &state) {
    for (int idx = 0; idx < (1 << state.n); ++idx) {
        if (abs(state.vals[idx]) < 1e-6) {
            state.vals[idx] = 0;
        }
        else {
            state.vals[idx] /= sqrt(abs(state.vals[idx]));
            state.vals[idx] = to_root8(from_root8(state.vals[idx]));
        }
    }
}

BruteState ground_brute(int n) {
    BruteState res;
    res.n = n;
    res.vals = std::vector<Complex>(1 << n, 0);
    res.vals[0] = 1;
    return res;
}

void brute_x(BruteState &state, int i) {
    for (int idx = 0; idx < (1 << state.n); ++idx)
        if (idx & (1 << i))
            std::swap(state.vals[idx], state.vals[idx ^ (1 << i)]);
    clean(state);
}

void brute_z(BruteState &state, int i) {
    for (int idx = 0; idx < (1 << state.n); ++idx)
        if (idx & (1 << i))
            state.vals[idx]*= -1;
    clean(state);
}

void brute_s(BruteState &state, int i) {
    for (int idx = 0; idx < (1 << state.n); ++idx)
        if (idx & (1 << i))
            state.vals[idx] *= I;
    clean(state);
}

void brute_cz(BruteState &state, int i, int j) {
    for (int idx = 0; idx < (1 << state.n); ++idx)
        if ((idx & (1 << i)) && (idx & (1 << j)))
            state.vals[idx] *= -1;
    clean(state);
}

void brute_swap(BruteState &state, int i, int j) {
    for (int idx = 0; idx < (1 << state.n); ++idx) {
        if ((idx & (1 << i)) && (~idx & (1 << j)))
            std::swap(state.vals[idx], state.vals[idx ^ (1 << i) ^ (1 << j)]);
    }
}

void brute_h(BruteState &state, int i) {
    std::vector<Complex> new_vals(1 << state.n);
    for (int idx = 0; idx < (1 << state.n); ++idx) {
        if (~idx & (1 << i)) {
            new_vals[idx & (~ (1 << i))]+= state.vals[idx] / sqrt(2);
            new_vals[idx | (1 << i)]+= state.vals[idx] / sqrt(2);
        }
        else {
            new_vals[idx & (~ (1 << i))]+= state.vals[idx] / sqrt(2);
            new_vals[idx | (1 << i)]-= state.vals[idx] / sqrt(2);
        }
    }

    state.vals = new_vals;
    clean(state);
}

void brute_cx(BruteState &state, int i, int j) {
    brute_h(state, j);
    brute_cz(state, i, j);
    brute_h(state, j);
}

void brute_push(BruteState &state) {
    state.n++;
    state.vals.resize(1 << state.n, 0);
    clean(state);
}

void brute_pop(BruteState &state) {
    std::vector<Complex> new_vals(1 << state.n - 1, 0);
    for (int idx = 0; idx < (1 << state.n); idx++) {
        if (~idx & (1 << state.n - 1))
            new_vals[idx & (~ (1 << state.n - 1))] = state.vals[idx];
    }
    state.vals = new_vals;
    state.n--;
    clean(state);
}

bool brute_eq(BruteState &state1, BruteState &state2) {
    if (state1.n != state2.n)
        return false;
    for (int idx = 0; idx < (1 << state1.n); ++idx)
        if (abs(state1.vals[idx] - state2.vals[idx]) > 1e-6)
            return false;
    return true;
}

bool brute_is_zero(BruteState &state) {
    for (int idx = 0; idx < (1 << state.n); ++idx)
        if (abs(state.vals[idx]) > 1e-6)
            return false;
    return true;
}

bool brute_proj_eq(BruteState &state1, BruteState &state2) {
    if (state1.n != state2.n)
        return false;

    bool difference_found = false;
    int phasediff = 0;

    for (int idx = 0; idx < (1 << state1.n); ++idx) {
        if (abs(state1.vals[idx]) < 1e-6 && abs(state2.vals[idx]) < 1e-6)
            continue;
        if (abs(state1.vals[idx]) < 1e-6 || abs(state2.vals[idx]) < 1e-6)
            return false;
        
        if (!difference_found) {
            difference_found = true;
            phasediff = from_root8(state1.vals[idx] / state2.vals[idx]);
        }

        if (from_root8(state1.vals[idx] / state2.vals[idx]) != phasediff)
            return false;
    }

    return true;
}

bool test_proj_eq(StabState state1, BruteState state2) {
    if (state1.n != state2.n)
        return false;

    if (state1.is_zero || brute_is_zero(state2))
        return state1.is_zero == brute_is_zero(state2);

    int N = state1.n;

    bool found_difference = false;
    int phasediff = 0;

    clean(state2);

    std::function<int(int)> get_stab_phase = [&](int idx) -> int {
        int phase = 0;
        for (int i = 0; i < N; ++i)
            for (int j = i + 1; j < N; ++j)
                if ((idx & (1 << i)) && (idx & (1 << j)))
                    phase+= 2 * int(state1.quad_part.get(i, j));

        for (int i = 0; i < N; ++i)
            if (idx & (1 << i))
                phase+= state1.lin_part[i];

        return 2 * (phase % 4);
    };

    for (int idx = 0; idx < (1 << N); ++idx) {
        bool in_support1, in_support2;

        BVector vec_idx(N);
        vec_idx.vec[0] = idx;

        in_support1 = state1.A * vec_idx == state1.b;
        in_support2 = abs(state2.vals[idx]) > 1e-6;

        if (in_support1 ^ in_support2)
            return false;
        if (!in_support1 && !in_support2)
            continue;

        if (!found_difference) {
            found_difference = true;
            phasediff = from_root8(to_root8(get_stab_phase(idx)) / state2.vals[idx]);
        }

        if (from_root8(to_root8(get_stab_phase(idx)) / state2.vals[idx]) != phasediff)
            return false;
    }

    return true;
}

bool test_eq(StabState state1, BruteState state2) {
    if (state1.n != state2.n)
        return false;

    if (state1.is_zero || brute_is_zero(state2))
        return state1.is_zero == brute_is_zero(state2);

    int N = state1.n;

    bool found_difference = false;
    int phasediff = 0;

    clean(state2);

    std::function<int(int)> get_stab_phase = [&](int idx) -> int {
        int phase = 0;
        for (int i = 0; i < N; ++i)
            for (int j = i + 1; j < N; ++j)
                if ((idx & (1 << i)) && (idx & (1 << j)))
                    phase+= 2 * int(state1.quad_part.get(i, j));

        for (int i = 0; i < N; ++i)
            if (idx & (1 << i))
                phase+= state1.lin_part[i];

        return 2 * (phase % 4);
    };

    for (int idx = 0; idx < (1 << N); ++idx) {
        bool in_support1, in_support2;

        BVector vec_idx(N);
        vec_idx.vec[0] = idx;

        in_support1 = state1.A * vec_idx == state1.b;
        in_support2 = abs(state2.vals[idx]) > 1e-6;

        if (in_support1 ^ in_support2)
            return false;
        if (!in_support1 && !in_support2)
            continue;

        if (!found_difference) {
            found_difference = true;
            phasediff = from_root8(to_root8(get_stab_phase(idx)) / state2.vals[idx]);
        }

        if (from_root8(to_root8(get_stab_phase(idx)) / state2.vals[idx]) != phasediff)
            return false;
    }

    return phasediff == (8 - state1.phase) % 8;
}

void print_brute(BruteState state) {
    Complex first_nonzero_val = 0;
    clean(state);
    for (int idx = 0; idx < (1 << state.n); ++idx) {
        if (abs(state.vals[idx]) > 1e-6) {
            first_nonzero_val = state.vals[idx];
            break;
        }
    }
    
    for (int idxt = 0; idxt < (1 << state.n); ++idxt) {
        int idx = 0;
        for (int i = 0; i < state.n; ++i) {
            if (idxt & (1 << i)) {
                idx |= (1 << (state.n - 1 - i));
            }
        }

        if (abs(state.vals[idx]) < 1e-6)
            continue;
        int phase = from_root8(state.vals[idx] / first_nonzero_val) / 2;
        bool wtf = from_root8(state.vals[idx] / first_nonzero_val) % 2;
        
        my_assert(!wtf);

        switch (phase) {
            case 0:
                std::cout << "  +";
                break;
            case 1:
                std::cout << " +i";
                break;
            case 2:
                std::cout << "  -";
                break;
            case 3:
                std::cout << " -i";
                break;
        }
        std::cout << "|";
        for (int i = 0; i < state.n; ++i)
            std::cout << ((idx & (1 << i)) ? 1 : 0);
        std::cout << ">";
    }
    std::cout << std::endl;
}

bool brute_test(int seed, int N, int depth0) {
    srand(seed);

    std::vector<Gate> gates = rand_circ(N, depth0);
    std::vector<BruteState> brute_states;
    std::vector<StabState> stab_states;

    brute_states.push_back(ground_brute(N));
    stab_states.push_back(ground_state(N));

    for (int ptr = 0; ptr < gates.size(); ++ptr) {
        auto gate = gates[ptr];

        brute_states.push_back(brute_states[ptr]);
        stab_states.push_back(stab_states[ptr]);

        if (gate.name == "H") {
            brute_h(brute_states[ptr + 1], gate.i);
            apply_h(stab_states[ptr + 1], gate.i);
        }
        else if (gate.name == "CZ") {
            brute_cz(brute_states[ptr + 1], gate.i, gate.j);
            apply_cz(stab_states[ptr + 1], gate.i, gate.j);
        }
        else if (gate.name == "SWAP") {
            brute_swap(brute_states[ptr + 1], gate.i, gate.j);
            apply_swap(stab_states[ptr + 1], gate.i, gate.j);
        }
        else if (gate.name == "X") {
            brute_x(brute_states[ptr + 1], gate.i);
            apply_x(stab_states[ptr + 1], gate.i);
        }        
        else if (gate.name == "Z") {
            brute_z(brute_states[ptr + 1], gate.i);
            apply_z(stab_states[ptr + 1], gate.i);
        }
        else if (gate.name == "S") {
            brute_s(brute_states[ptr + 1], gate.i);
            apply_s(stab_states[ptr + 1], gate.i);
        }
        else if (gate.name == "PUSH") {
            push_qubit(stab_states[ptr + 1]);
            brute_push(brute_states[ptr + 1]);
        }
        else if (gate.name == "POP") {
            pop_qubit(stab_states[ptr + 1]);
            brute_pop(brute_states[ptr + 1]);
        }
        
        bool brute_zero = brute_is_zero(brute_states[ptr + 1]);
        bool stab_zero = stab_states[ptr + 1].is_zero;

        if (brute_zero || stab_zero)
            break;

        stab_states[ptr + 1] = normal_form(stab_states[ptr + 1]);        
    }


    for (int ptr = 1; ptr < brute_states.size(); ++ptr) {
        if (!test_eq(stab_states[ptr], brute_states[ptr])) {
            std::cout << "Failed with seed " << seed << ", depth " << depth0 << " and " << N << " qubits" << std::endl;
            print_brute(brute_states[ptr - 1]);
            std::cout << "Applied " << gates[ptr - 1].name << " " << gates[ptr - 1].i;
            if (gates[ptr - 1].name == "S" || gates[ptr - 1].name == "H" || gates[ptr - 1].name == "X" || gates[ptr - 1].name == "Z")
                std::cout << std::endl;
            else
                std::cout << " " << gates[ptr - 1].j << std::endl;
            std::cout << "Brute result: " << std::endl;
            print_brute(brute_states[ptr]);
            std::cout << "Stab result: " << std::endl;
            print_superposition(stab_states[ptr]);
            print(stab_states[ptr]);
            std::cout << "--------------------------------" << std::endl;
            return true;
        }
    }

    return false;
}

bool brute_test_suite() {
    std::vector<int> depths;
    depths.push_back(10);
    depths.push_back(20);
    depths.push_back(50);
    depths.push_back(100);

    for (int N = 3; N <= 10; ++N) {
        for (int depth : depths) {
            std::cout << "Testing " << N << " qubits with depth " << depth << std::endl;
            for (int it = 0; it < 1000; ++it) {
                if (brute_test(it, N, depth))
                    return true;
            }
            std::cout << "Success!" << std::endl;
        }
    }
    return false;
}



int main() {
    brute_test_suite();
//    std::cout << clifford_group_size(2) << std::endl;

    return 0;
    auto res = ground_state(4);
    apply_h(res, 0);
    apply_cx(res, 0, 2);
    apply_h(res, 2);

    apply_h(res, 1);
    apply_cx(res, 1, 3);
    apply_s(res, 1);
    res = normal_form(res);

    apply_cx(res, 2, 3);
    res = normal_form(res);

    apply_h(res, 2); // fault is here
    res = normal_form(res);

    print(res);
    print_superposition(res);

    pop_qubit(res);
    pop_qubit(res);
    res = normal_form(res);
    print(res);
    print_superposition(res);
    std::cout << "--------------------------------" << std::endl;

    return 0;
}
