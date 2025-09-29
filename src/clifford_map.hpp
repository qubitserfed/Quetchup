#pragma once

#include "clifford_state.hpp"

#include <algorithm>
#include <vector>

struct CliffordMap {
    StabState state;
    int in_wires, out_wires;
};

bool operator == (CliffordMap f, CliffordMap g);
CliffordMap compose(CliffordMap f, CliffordMap g);
CliffordMap tensor(CliffordMap f, CliffordMap g);
CliffordMap id_map(int n);
CliffordMap h_map();
CliffordMap x_map();
CliffordMap y_map();
CliffordMap z_map();
CliffordMap i_map();
CliffordMap j_map();
CliffordMap s_map();
CliffordMap cx_map();
CliffordMap cz_map();
CliffordMap zero_projector();
CliffordMap swap_map();
CliffordMap from_state(StabState state);
StabState   apply_map(CliffordMap f, StabState state);
