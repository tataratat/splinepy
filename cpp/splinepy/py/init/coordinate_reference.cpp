#include <splinepy/py/py_coordinate_references.hpp>

namespace splinepy::py::init {
namespace py = pybind11;
void init_coordinate_references(py::module_& m) {
  add_coordinate_references_pyclass(m);
}

} // namespace splinepy::py::init
