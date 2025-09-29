#pragma once

#include <algorithm>
#include <vector>
#include <cassert>

using u64 = unsigned long long;

const u64 x_mask = 0xAAAAAAAAAAAAAAAAULL;
const u64 z_mask = ~x_mask;

struct BVector {
    std::vector<u64> vec;
    int n;

    int  no_buckets();
    bool get(int);
    void flip(int);
    bool is_zero();
    void set(int , bool);
    int  weight();
    void push_back(bool);
    void pop_back();
    void swap(int, int);

    BVector();
    BVector(int);
    BVector(std::vector<bool>);
};


struct BMatrix {
    std::vector<BVector> mat;
    int n, m;

    BVector &row(int);
    BVector column(int);
    BVector &last_row();

    void append_row(BVector);
    void flip(int, int);
    void pop_row();
    void set(int, int, bool);
    bool get(int, int);
    bool empty();

    void append_column(BVector);
    void pop_column();
    void add_rows(int, int);
    void swap_rows(int, int);
    void swap_cols(int, int);
    void remove_zeros();
    void sort_rows();

    BMatrix row_submatrix(std::vector<int>);
    BMatrix column_submatrix(std::vector<int>);

    BMatrix();
    BMatrix(int, int);
};

void                my_assert                       (bool arg);
int                 popcount                        (u64 num);

// BVector definitions
void                print                           (BVector vec);
BVector             operator +                      (const BVector &, const BVector &);
bool                operator *                      (const BVector &, const BVector &);
bool                operator <                      (const BVector &, const BVector &);
bool                operator ==                     (const BVector &, const BVector &);
bool                operator !=                     (const BVector &, const BVector &);
bool                sym_prod                        (BVector, BVector);
BVector             versor                          (int n, int i);

// BMatrix definitions
void                print                           (BMatrix);
BMatrix             identity                        (const int &);
BMatrix             transpose                       (BMatrix);
BVector             operator *                      (BVector, BMatrix);
BVector             operator *                      (BMatrix, BVector);
BMatrix             operator *                      (BMatrix, BMatrix);
bool                operator ==                     (BMatrix, BMatrix);
bool                operator !=                     (BMatrix, BMatrix);
bool                is_zero                         (BMatrix);

// Linear algebra algorithms
int                 to_row_echelon                  (BMatrix &);
std::vector<int>    restricted_row_echelon          (BMatrix &, std::vector<int>);
void                to_row_echelon                  (BMatrix &, BVector&);
std::vector<int>    get_pivots                      (BMatrix &);
std::vector<int>    get_pivots                      (BMatrix &, std::vector<int>);
BVector             canonical_quotient              (BVector, BMatrix &);
bool                in_span                         (BMatrix, const BVector &);
BMatrix             transposed_product              (BMatrix &, BMatrix &);
BVector             transposed_product              (const BVector &, BMatrix &);
BMatrix             basis_completion                (BMatrix);
constexpr u64       raw_sym_prod_ll                 (const u64 &, const u64&);
BMatrix             isotropic_closure               (BMatrix);
BMatrix             find_kernel                     (BMatrix);
BVector             solve                           (BMatrix, BVector);

std::pair<BMatrix,BVector> affine_extension        (BMatrix, BVector, BVector);