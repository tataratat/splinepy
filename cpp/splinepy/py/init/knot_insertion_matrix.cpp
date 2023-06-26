#include <splinepy/py/py_knot_insertion_matrix.hpp>

namespace splinepy::py::init {

void init_knot_insertion_matrix(py::module_& m) {
  add_knot_insertion_matrix(m);
}

} // namespace splinepy::py::init
