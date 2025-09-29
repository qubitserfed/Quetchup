#pragma once

#include <string>
#include "linear_algebra.hpp"

struct StabState {
    int n;

    int phase; // global phase, not used for now
    int magnitude; // magnitude of the state, not used for now
    bool is_zero; // whether the state is the zero state

    std::vector<int> lin_part;
    BMatrix quad_part; // diagonal gets ignored

    StabState(int n) {
        this->n = n;
        this->phase = 0;
        this->magnitude = 0;
        this->is_zero = false;
        this->lin_part = std::vector<int>(n, 0);
        this->quad_part = BMatrix(n, n);

        this->A = identity(n);
        this->b = BVector(n);
    }

    StabState() {
        StabState(0);
    }

    BMatrix A;
    BVector b; // define affine space Ax=b
};

bool operator == (StabState, StabState);

StabState   ground_state        (int);
StabState   bell_states         (int);
StabState   normal_form         (StabState);
void        apply_cz            (StabState &, int, int);
void        apply_cx            (StabState &, int, int);
void        apply_swap          (StabState &, int, int);
void        apply_x             (StabState &, int);
void        apply_z             (StabState &, int);
void        apply_h             (StabState &, int);
void        apply_s             (StabState &, int);
void        push_qubit          (StabState &);
void        pop_qubit           (StabState &);
StabState   tensor              (StabState, StabState);
void        permute             (StabState &, std::vector<int>);
void        print               (StabState);
void        print_compact       (StabState);
bool        is_ground           (StabState);
void        print_superposition (StabState);
bool        real_proj_eq        (StabState, StabState);
std::string to_string           (StabState);
std::string to_latex            (StabState);
