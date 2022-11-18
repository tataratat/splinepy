#include <splinepy/py/py_spline_exporter.hpp>

namespace splinepy::py::init {

void init_exporter(py::module_& m) {
  // Void functions that define arguments
  // returns [connectivity, vertex_ids, edge_information, boundaries]
  m.def("retrieve_mfem_information",
        &splinepy::py::RetrieveMfemInformation,
        py::arg("corner_vertices"),
        py::arg("tolerance"));
  m.def("export_iges",
        &splinepy::py::ExportIges,
        py::arg("fname"),
        py::arg("splines"));
  m.def("export_irit",
        &splinepy::py::ExportIrit,
        py::arg("fname"),
        py::arg("splines"));
  m.def("export_xml",
        &splinepy::py::ExportXml,
        py::arg("fname"),
        py::arg("splines"));
  m.def("export_vtk",
        &splinepy::py::ExportVtk,
        py::arg("fname"),
        py::arg("splines"),
        py::arg("resolutions_per_spline"));
}

} // namespace splinepy::py::init
