#include <splinepy/py/py_spline_reader.hpp>

namespace splinepy::py::init {

void init_reader(py::module_& m) { AddSplineReader(m); }

} // namespace splinepy::py::init
