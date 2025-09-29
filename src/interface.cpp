#include <iostream>

#include "clifford_state.hpp"

int main(int argc, char **argv) {
    std::string gate;
    bool verbose_flag;
    int src, targ;
    int n;


    verbose_flag = false;
    if (argc > 1) {
        for (int i = 1; i < argc; ++i) {
            if (std::string(argv[i]) == "-v") {
                verbose_flag = true;
            }
        }
    }

    std::cin >> n;
    StabState state = ground_state(n);

    while (std::cin >> gate) {
        if (gate == "H") {
            std::cin >> targ;
            apply_h(state, targ);
        }
        else if (gate == "X") {
            std::cin >> targ;
            apply_x(state, targ);
        }
        else if (gate == "Z") {
            std::cin >> targ;
            apply_z(state, targ);
        }
        else if (gate == "S") {
            std::cin >> targ;
            apply_s(state, targ);
        }
        else if (gate == "CNOT") {
            std::cin >> src >> targ;
            apply_cx(state, src, targ);
        }
        else if (gate == "SWAP") {
            std::cin >> src >> targ;
            apply_swap(state, src, targ);
        }
        else if (gate == "CZ") {
            std::cin >> src >> targ;
            apply_cz(state, src, targ);
        }
        else if (gate == "PUSH") {
            push_qubit(state);
        }
        else if (gate == "POP") {
            pop_qubit(state);
        }
        state = normal_form(state);
    }

    if (verbose_flag) {
        print(state);
        print_superposition(state);
    }
    else
        print_compact(state);


    return 0;
}