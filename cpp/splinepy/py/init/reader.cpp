#include <splinepy/py/spline_reader.hpp>

void init_reader(py::module_& m) {
  // Functions that return list of dict.
  // Keys are ["knot_vectors", "control_points", "degrees"] (+ ["weights"])
  m.def("read_iges", &read_iges, py::arg("fname"))
      .def("read_xml", &read_xml, py::arg("fname"))
      .def("read_irit", &read_irit, py::arg("fname"));
}
