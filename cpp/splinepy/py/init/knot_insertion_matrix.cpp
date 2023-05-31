#include <splinepy/py/py_knot_insertion_matrix.hpp>

namespace splinepy::py::init {

void init_knot_insertion_matrix(py::module_& m) {
  // Functions that return fitted bspline as dict.
  m.def("knot_insertion_matrix",
        &splinepy::py::ComputeKnotInsertionMatrix,
        py::arg("old_knot_vector"),
        py::arg("new_knot_vector"),
        py::arg("degree"),
        py::arg("tolerance"));
  m.def("global_knot_insertion_matrix",
        &splinepy::py::ComputeGlobalKnotInsertionMatrix,
        py::arg("old_knot_vectors"),
        py::arg("new_knot_vectors"),
        py::arg("degrees"),
        py::arg("tolerance"));
}

} // namespace splinepy::py::init
