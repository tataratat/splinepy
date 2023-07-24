#include "splinepy/py/py_knot_vector.hpp"

namespace splinepy::py::init {

void init_knot_vector(py::module_& m) { add_knot_vector(m); }

} // namespace splinepy::py::init
