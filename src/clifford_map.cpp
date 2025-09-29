#include "clifford_map.hpp"
#include "clifford_state.hpp"
#include "linear_algebra.hpp"

#include <iostream>


bool operator == (CliffordMap f, CliffordMap g) {
    return f.state == g.state && f.in_wires == g.in_wires && f.out_wires == g.out_wires;
}

CliffordMap tensor(CliffordMap f, CliffordMap g) {
    CliffordMap h;

    h.state = tensor(f.state, g.state);
    h.in_wires = f.in_wires + g.in_wires;
    h.out_wires = f.out_wires + g.out_wires;

    std::vector<int> perm(h.state.n);
    // product in_wires
    for (int i = 0; i < f.in_wires; ++i)
        perm[i] = i;
    for (int i = 0; i < g.in_wires; ++i)
        perm[i + f.in_wires] = i + f.in_wires + f.out_wires;
    
    for (int i = 0; i < f.out_wires; ++i)
        perm[i + f.in_wires + g.in_wires] = i + f.in_wires;
    for (int i = 0; i < g.out_wires; ++i)
        perm[i + f.in_wires + g.in_wires + f.out_wires] = i + f.in_wires + f.out_wires + g.in_wires;

    permute(h.state, perm);
    return h;
}

CliffordMap compose(CliffordMap f, CliffordMap g) {
    if (f.out_wires != g.in_wires) {
        my_assert(0);
    }


    StabState h_state = tensor(f.state, g.state);

    std::vector<int> perm(h_state.n);
    for (int i = 0; i < f.in_wires; ++i)
        perm[i] = i;

    for (int i = 0; i < g.out_wires; ++i)
        perm[i + f.in_wires] = i + f.in_wires + f.out_wires + g.in_wires;

    for (int i = 0; i < g.in_wires; ++i) {
        perm[2 * i + f.in_wires + g.out_wires] = f.in_wires + i;
        perm[2 * i + 1 + f.in_wires + g.out_wires] = f.in_wires + f.out_wires + i;
    }

    std::vector<int> inv_perm(h_state.n);
    for (int i = 0; i < h_state.n; ++i)
        inv_perm[perm[i]] = i;
    perm = inv_perm;


    permute(h_state, perm); h_state = normal_form(h_state);

    for (int i = 0; i < f.out_wires; ++i) {
        apply_cx(h_state, h_state.n - 2 * i - 2, h_state.n - 2 * i - 1);
        h_state = normal_form(h_state);
        apply_h(h_state, h_state.n - 2 * i - 2);
        h_state = normal_form(h_state);
    }

    for (int i = 0; i < f.out_wires; ++i) {
        pop_qubit(h_state);
        pop_qubit(h_state);
        h_state.magnitude-= 1;
    }

    CliffordMap h;
    h.state = normal_form(h_state);
    h.in_wires = f.in_wires;
    h.out_wires = g.out_wires;
    return h;
}

CliffordMap id_map(int n) {
    CliffordMap h;
    h.state = ground_state(2 * n);
    h.in_wires = n;
    h.out_wires = n;
    h.state.magnitude = -n;
    for (int i = 0; i < n; ++i) {
        apply_h(h.state, i);
        apply_cx(h.state, i, i + n);
    }
    h.state = normal_form(h.state);
    
    return h;
}

CliffordMap swap_map() {
    CliffordMap h = id_map(2);
    apply_swap(h.state, 2, 3);
    h.state = normal_form(h.state);
    return h;
}

CliffordMap h_map() {
    CliffordMap h = id_map(1);
    apply_h(h.state, 1);
    h.state = normal_form(h.state);
    return h;
}

CliffordMap x_map() {
    CliffordMap h = id_map(1);
    apply_x(h.state, 1);
    h.state = normal_form(h.state);
    return h;
}

CliffordMap z_map() {
    CliffordMap h = id_map(1);
    apply_z(h.state, 1);
    h.state = normal_form(h.state);
    return h;
}

CliffordMap y_map() {
    CliffordMap h = id_map(1);
    apply_z(h.state, 1);
    apply_x(h.state, 1);
    h.state.phase = (h.state.phase + 2) % 8;
    return h;
}

CliffordMap i_map() {
    CliffordMap h = id_map(0);
    h.state.phase = 2;
    h.in_wires = 0, h.out_wires = 0;
    return h;
}

CliffordMap j_map() {
    CliffordMap h = id_map(0);
    h.state.phase = 1;
    h.in_wires = 0, h.out_wires = 0;
    return h;
}

CliffordMap s_map() {
    CliffordMap h = id_map(1);
    apply_s(h.state, 1);
    h.state = normal_form(h.state);
    return h;
}

CliffordMap cx_map() {
    CliffordMap h = id_map(2);
    apply_cx(h.state, 2, 3);
    h.state = normal_form(h.state);
    return h;
}

CliffordMap cz_map() {
    CliffordMap h = id_map(2);
    apply_cz(h.state, 2, 3);
    h.state = normal_form(h.state);
    return h;
}

CliffordMap zero_projector() {
    CliffordMap h = id_map(1);
    push_qubit(h.state);
    apply_cx(h.state, 1, 2);
    pop_qubit(h.state);
    h.state.magnitude+= 1;
    h.in_wires = 1, h.out_wires = 0;
    return h;
}

CliffordMap from_state(StabState state) {
    CliffordMap h;
    h.state = state;
    h.in_wires = 0;
    h.out_wires = state.n;
    return h;
}

StabState apply_map(CliffordMap f, StabState state) {
    if (f.in_wires != state.n) {
        std::cerr << "The number of input wires of the map and the size of the state do not match" << std::endl;
        my_assert(0);
    }

    StabState res = tensor(state, f.state);

    for (int i = 0; i < state.n; ++i) {
        apply_cx(res, i, i + state.n);
        apply_h(res, i);
    }
    
    std::vector<int> perm(res.n);

    // place state wires
    for (int i = 0; i < state.n; ++i) // [0, in_wires) -> [2 * in_wires, 2 * in_wires + state.n)
        perm[i] = i + 2 * state.n;

    // place in wires
    for (int i = 0; i < f.in_wires; ++i) // [in_wires, in_wires + f.in_wires) -> [f.in_wires, 2 * f.in_wires)
        perm[i + state.n] = i + state.n;

    //place out wires
    for (int i = 0; i < f.out_wires; ++i) // [2 * in_wires, 2 * in_wires + f.out_wires) -> [state.n, state.n + f.out_wires)
        perm[i + 2 * state.n] = i;

    std::vector<int> inv_perm(res.n);
    for (int i = 0; i < res.n; ++i)
        inv_perm[perm[i]] = i;
    perm = inv_perm;

    permute(res, perm);
    for (int i = 0; i < 2 * state.n; ++i)
        pop_qubit(res);

    res.magnitude-= f.in_wires;
    return res;
}
