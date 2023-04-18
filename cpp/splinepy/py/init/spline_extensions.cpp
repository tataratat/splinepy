#include <splinepy/py/py_spline_extensions.hpp>

namespace splinepy::py::init {
namespace py = pybind11;
void init_spline_extensions(py::module_& m) { add_spline_extensions(m); }

} // namespace splinepy::py::init
