#include <splinepy/py/spline_exporter.hpp>

void init_exporter(py::module_& m) {
  // Void functions that define arguments
  // returns [connectivity, vertex_ids, edge_information, boundaries]
  m.def("retrieve_mfem_information",
        &retrieve_MFEM_information,
        py::arg("corner_vertices"),
        py::arg("tolerance"));
}
