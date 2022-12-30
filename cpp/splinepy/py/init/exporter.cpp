#include <splinepy/py/py_spline_exporter.hpp>

namespace splinepy::py::init {

void init_exporter(py::module_& m) {
  // Void functions that define arguments
  // returns [connectivity, vertex_ids, edge_information, boundaries]
  m.def("retrieve_mfem_information",
        &splinepy::py::RetrieveMfemInformation,
        py::arg("corner_vertices"),
        py::arg("tolerance"));
  m.def("interfaces_from_boundary_centers",
        &splinepy::py::InterfacesFromBoundaryCenters,
        py::arg("face_center_vertices"),
        py::arg("tolerance"),
        py::arg("para_dim"));
  m.def("get_orientation",
        &splinepy::py::GetBoundaryOrientation,
        py::arg("start_spline"),
        py::arg("start_boundary_id"),
        py::arg("end_spline"),
        py::arg("end_boundary_id"),
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
