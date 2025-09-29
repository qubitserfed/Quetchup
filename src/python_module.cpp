#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "clifford_state.hpp"
#include "clifford_map.hpp"

namespace py = pybind11;

PYBIND11_MODULE(Extension, module) {
    module.doc() = "Qubos is a library for quantum Stabilizer simulation";

    py::class_<CliffordMap>(module, "Map")
        .def(py::init<>())
        .def("__eq__", static_cast<bool (*)(CliffordMap, CliffordMap)>(&operator==))
        .def("__mul__", &compose)
        .def("magnitude", [&](CliffordMap &self) -> int {
            return self.state.magnitude;
        })
        .def("phase", [&](CliffordMap &self) -> int {
            return self.state.phase;
        })
        .def("in_wires", [&](CliffordMap &self) -> int {
            return self.in_wires;
        })
        .def("out_wires", [&](CliffordMap &self) -> int {
            return self.out_wires;
        })
        .def("s", [&](CliffordMap& self, int pos) -> CliffordMap& {
            apply_s(self.state, pos + self.in_wires);
            return self;
        })
        .def("h", [&](CliffordMap& self, int pos) -> CliffordMap& {
            apply_h(self.state, pos + self.in_wires);
            return self;
        })
        .def("x", [&](CliffordMap& self, int pos) -> CliffordMap& {
            apply_x(self.state, pos + self.in_wires);
            return self;
        })
        .def("z", [&](CliffordMap& self, int pos) -> CliffordMap& {
            apply_z(self.state, pos + self.in_wires);
            return self;
        })
        .def("cx", [&](CliffordMap& self, int pos0, int pos1) -> CliffordMap& {
            apply_cx(self.state, pos0 + self.in_wires, pos1 + self.in_wires);
            return self;
        })
        .def("cz", [&](CliffordMap& self, int pos0, int pos1) -> CliffordMap& {
            apply_cz(self.state, pos0 + self.in_wires, pos1 + self.in_wires);
            return self;
        })
        .def("swap", [&](CliffordMap& self, int pos0, int pos1) -> CliffordMap& {
            apply_swap(self.state, pos0 + self.in_wires, pos1 + self.in_wires);
            return self;
        })
        .def("prep_zero", [&](CliffordMap& self) -> CliffordMap& {
            push_qubit(self.state);
            self.out_wires+= 1;
            return self;
        })
        .def("prep_zero", [&](CliffordMap& self, int pos) -> CliffordMap& {
            push_qubit(self.state);
            for (int i = self.state.n - 1; i > pos + self.in_wires; --i)
                apply_swap(self.state, i, i - 1);
            self.out_wires+= 1;
            return self;
        })
        .def("post_zero", [&](CliffordMap& self, int pos) -> CliffordMap& {
            for (int i = pos + self.in_wires; i < self.state.n - 1; ++i)
                apply_swap(self.state, i, i + 1);
            pop_qubit(self.state);
            self.out_wires-= 1;
            return self;
        })
        .def("clone", [&](CliffordMap self) -> CliffordMap {
            return self;
        });

    module.def("map_tensor", static_cast<CliffordMap (*)(CliffordMap, CliffordMap)>(&tensor));
    module.def("compose", &compose);
    module.def("id_map", &id_map);
    module.def("h_map", &h_map); 
    module.def("x_map", &x_map);
    module.def("z_map", &z_map);
    module.def("y_map", &y_map);
    module.def("i_map", &i_map);
    module.def("j_map", &j_map);
    module.def("s_map", &s_map);
    module.def("cx_map", &cx_map);
    module.def("cz_map", &cz_map);
    module.def("swap_map", &swap_map);
    module.def("zero_post", &zero_projector);
    module.def("zero_prep", [&]() -> CliffordMap {
        CliffordMap res;
        push_qubit(res.state);
        res.in_wires = 0, res.out_wires = 1;
        return res;
    });

    py::class_<StabState>(module, "State")
        .def(py::init<>())
        .def(py::init<int>())
        .def("__eq__", static_cast<bool (*)(StabState, StabState)>(&operator==))
        .def("magnitude", [&](StabState self) -> int {
            normal_form(self);
            return self.magnitude + self.n - self.A.n;
        })
        .def("phase", [&](StabState self) -> int { return self.phase; })
        .def("n", [&](StabState self) -> int { return self.n; })
        .def("is_zero", [&](StabState self) -> bool { return self.is_zero; })
        .def("affine_part", [&](StabState self) -> py::tuple {
            py::array_t<int> A = py::array_t<int>({self.A.n, self.A.m});
            auto buffA = A.mutable_unchecked<2>();
            for (int i = 0; i < self.A.n; ++i) {
                for (int j = 0; j < self.A.m; ++j) {
                    buffA(i, j) = self.A.get(i, j);
                }
            }

            py::array_t<int> b(self.b.n);
            auto buffb = b.mutable_unchecked<1>();
            for (int i = 0; i < self.b.n; ++i) {
                buffb(i) = self.b.get(i) ? 1 : 0;
            }

            return py::make_tuple(A, b);
        })
        .def("phase_polynomial_matrix", [&](StabState self) -> py::array_t<int> {
            py::array_t<int> res = py::array_t<int>({self.n, self.n});
            auto buff = res.mutable_unchecked<2>();
            for (int i = 0; i < self.n; ++i) {
                for (int j = 0; j < self.n; ++j) {
                    buff(i, j) = self.quad_part.get(i, j);
                }
            }
            return res;
        })
        .def("__latex_rep", [&](StabState self) -> std::string {
            return to_latex(self);
        })
        .def("prep_zero", [&](StabState &self) -> StabState {
            push_qubit(self);
            return self;
        })
        .def("prep_zero", [&](StabState &self, int pos) -> StabState {
            push_qubit(self);
            for (int i = self.n - 1; i > pos; --i)
                apply_swap(self, i, i - 1);
            return self;
        })
        .def("post_zero", [&](StabState &self, int pos) -> StabState {
            for (int i = pos; i < self.n - 1; ++i)
                apply_swap(self, i, i + 1);
            pop_qubit(self);
            return self;
        })
        .def("h", [&](StabState &self, int pos) -> StabState {
            apply_h(self, pos);
            return self;
        })
        .def("x", [&](StabState &self, int pos) -> StabState {
            apply_x(self, pos);
            return self;
        })
        .def("z", [&](StabState &self, int pos) -> StabState {
            apply_z(self, pos);
            return self;
        })
        .def("cx", [&](StabState &self, int pos0, int pos1) -> StabState {
            apply_cx(self, pos0, pos1);
            return self;
        })
        .def("cz", [&](StabState &self, int pos0, int pos1) -> StabState {
            apply_cz(self, pos0, pos1);
            return self;
        })
        .def("swap", [&](StabState &self, int pos0, int pos1) -> StabState {
            apply_swap(self, pos0, pos1);
            return self;
        })
        .def("s", [&](StabState &self, int pos) -> StabState {
            apply_s(self, pos);
            return self;
        })
        .def("clone", [&](StabState self) -> StabState {
            return self;
        });


    module.def("state_tensor", static_cast<StabState (*)(StabState, StabState)>(&tensor));
    module.def("ground_state", &ground_state);

    module.def("apply_cz", &apply_cz);
    module.def("apply_cx", &apply_cx);
    module.def("apply_swap", &apply_swap);
    module.def("apply_x", &apply_x);
    module.def("apply_z", &apply_z);
    module.def("apply_h", &apply_h);
    module.def("apply_s", &apply_s);
    module.def("push_qubit", &push_qubit);
    module.def("pop_qubit", &pop_qubit);
    module.def("apply_map", &apply_map);
    module.def("choi_state", [&](CliffordMap &self) -> StabState {
        return self.state;
    });

    module.def("to_map", &from_state);
}
