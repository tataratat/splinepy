#include <splinepy/py/py_spline_list.hpp>

namespace splinepy::py::init {
namespace py = pybind11;
void init_spline_list(py::module_& m) { add_spline_list_pyclass(m); }

} // namespace splinepy::py::init
