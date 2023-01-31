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
  m.def("extract_all_boundary_splines",
        &splinepy::py::ExtractAllBoundarySplines,
        py::arg("splines"),
        py::arg("interfaces"),
        py::arg("nthreads") = 1);
  m.def("orientations",
        &splinepy::py::GetBoundaryOrientations,
        py::arg("splines"),
        py::arg("base_ids"),
        py::arg("base_face_ids"),
        py::arg("neighbor_ids"),
        py::arg("neighbor_face_ids"),
        py::arg("tolerance"),
        py::arg("nthreads") = 1);
  m.def("boundaries_from_continuity",
        &splinepy::py::AddBoundariesFromContinuity,
        py::arg("boundary_splines"),
        py::arg("boundary_interfaces"),
        py::arg("global_interfaces"),
        py::arg("tolerance"),
        py::arg("nthreads") = 1);
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
