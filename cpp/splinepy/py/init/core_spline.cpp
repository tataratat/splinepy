#include <splinepy/py/py_spline.hpp>

namespace splinepy::py::init {

void init_core_spline(py::module_& m) {
  add_spline_pyclass(m, "CoreSpline");
}

} // namespace splinepy::py::init
