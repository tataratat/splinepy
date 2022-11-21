#include <splinepy/py/py_spline_reader.hpp>

namespace splinepy::py::init {

void init_reader(py::module_& m) {
  // Functions that return list of dict.
  // Keys are ["knot_vectors", "control_points", "degrees"] (+ ["weights"])
  m.def("read_iges", &read_iges, py::arg("fname"))
      .def("read_xml", &read_xml, py::arg("fname"))
      .def("read_irit", &read_irit, py::arg("fname"));
}

} // namespace splinepy::py::init
