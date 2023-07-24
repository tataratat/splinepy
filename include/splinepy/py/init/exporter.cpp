#include <splinepy/py/py_spline_exporter.hpp>

namespace splinepy::py::init {

void init_exporter(py::module_& m) { add_spline_exporter(m); }

} // namespace splinepy::py::init
