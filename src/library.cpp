// C ABI wrappers for exporting functions in clifford_state.hpp and clifford_map.hpp

#include <cstdlib>
#include <cstring>
#include <vector>

#include "clifford_state.hpp"
#include "clifford_map.hpp"

extern "C" {

// -------------------- StabState API --------------------

StabState *qb_ground_state(int n) {
    StabState res = ground_state(n);
    return new StabState(res);
}

StabState *qb_normal_form(const StabState *state) {
    StabState res = normal_form(*state);
    return new StabState(res);
}

void qb_apply_cz(StabState *state, int i, int j) {
    apply_cz(*state, i, j);
}

void qb_apply_cx(StabState *state, int i, int j) {
    apply_cx(*state, i, j);
}

void qb_apply_swap(StabState *state, int i, int j) {
    apply_swap(*state, i, j);
}

void qb_apply_x(StabState *state, int i) {
    apply_x(*state, i);
}

void qb_apply_z(StabState *state, int i) {
    apply_z(*state, i);
}

void qb_apply_h(StabState *state, int i) {
    apply_h(*state, i);
}

void qb_apply_s(StabState *state, int i) {
    apply_s(*state, i);
}

void qb_push_qubit(StabState *state) {
    push_qubit(*state);
}

void qb_pop_qubit(StabState *state) {
    pop_qubit(*state);
}

StabState *qb_tensor_state(StabState *u, StabState *v) {
    return new StabState(tensor(*u, *v));
}

void qb_permute(StabState *state, int *perm, int len) {
    std::vector<int> v(perm, perm + len);
    permute(*state, v);
}

void qb_is_ground(const StabState *state, bool *res) {
    *res = is_ground(*state);
}

void qb_real_proj_eq(StabState *u, StabState *v, bool *res) {
    *res = real_proj_eq(*u, *v);
}

void qb_stabstate_eq(StabState *u, StabState *v, bool *res) {
    *res = ((*u) == (*v));
}

void qb_free_stabstate(StabState *state) {
    delete state;
}

void qb_free_string(char *s) {
    std::free(s);
}

void qb_get_affine_space( StabState *state, bool **A, bool *b) {
    A = new bool*[state->A.n];
    for (int i = 0; i < state->A.n; ++i) {
        A[i] = new bool[state->A.m];
        for (int j = 0; j < state->A.m; ++j) {
            A[i][j] = state->A.get(i, j);
        }
    }

    b = new bool[state->b.n];
    for (int i = 0; i < state->b.n; ++i) {
        b[i] = state->b.get(i);
    }
}

void get_n(StabState *state, int *n) {
    *n = state->n;
}

void get_phase(StabState *state, int *phase) {
    *phase = state->phase;
}

void get_magnitude(StabState *state, int *magnitude) {
    *magnitude = state->magnitude;
}

void get_is_zero(StabState *state, bool *is_zero) {
    *is_zero = state->is_zero;
}

void get_phase_polynomial_matrix(StabState *state, int **phase_polynomial_matrix) {
    phase_polynomial_matrix = new int*[state->n];
    for (int i = 0; i < state->n; ++i) {
        phase_polynomial_matrix[i] = new int[state->n];
        std::fill(phase_polynomial_matrix[i], phase_polynomial_matrix[i] + state->n, 0);
    }
    for (int i = 0; i < state->n; ++i) {
        for (int j = 0; j < state->n; ++j) {
            phase_polynomial_matrix[i][j] = 2 * int(state->quad_part.get(i, j));
        }
    }
    for (int i = 0; i < state->n; ++i) {
        phase_polynomial_matrix[i][i] = state->lin_part[i];
    }
}

// -------------------- CliffordMap API --------------------
CliffordMap *qb_compose_map(CliffordMap *f, CliffordMap *g) {
    return new CliffordMap(compose(*f, *g));
}

CliffordMap *qb_tensor_map(CliffordMap *f, CliffordMap *g) {
    return new CliffordMap(tensor(*f, *g));
}

CliffordMap *qb_id_map(int n) {
    return new CliffordMap(id_map(n));
}

CliffordMap *qb_h_map() {
    return new CliffordMap(h_map());
}

CliffordMap *qb_x_map() {
    return new CliffordMap(x_map());
}

CliffordMap *qb_z_map() {
    return new CliffordMap(z_map());
}

CliffordMap *qb_s_map() {
    return new CliffordMap(s_map());
}

CliffordMap *qb_cx_map() {
    return new CliffordMap(cx_map());
}

CliffordMap *qb_cz_map() {
    return new CliffordMap(cz_map());
}

CliffordMap *qb_zero_prep() {
    return new CliffordMap { ground_state(1), 0, 1 };
}

CliffordMap *qb_zero_post() {
    return new CliffordMap { ground_state(1), 1, 0 };
}

void qb_free_clifford_map(CliffordMap *m) {
    delete m;
}

}
