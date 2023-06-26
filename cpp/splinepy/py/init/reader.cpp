#include <splinepy/py/py_spline_reader.hpp>

namespace splinepy::py::init {

void init_reader(py::module_& m) { add_spline_reader(m); }

} // namespace splinepy::py::init
