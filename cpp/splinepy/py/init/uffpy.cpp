#include <splinepy/py/py_uffpy.hpp>

namespace splinepy::py::init {

namespace py = pybind11;
void init_uffpy(py::module_& m) { add_uffpy(m); }

} // namespace splinepy::py::init
