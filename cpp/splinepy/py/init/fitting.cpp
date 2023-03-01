#include <splinepy/py/py_fitting.hpp>

namespace splinepy::py::init {

void init_fitting(py::module_& m) {
  // Functions that return fitted bspline as dict.
  m.def("interpolate_curve",
        &splinepy::py::InterpolateCurve,
        py::arg("points"),
        py::arg("degree"),
        py::arg("centripetal"),
        py::arg("knot_vector"));
  m.def("approximate_curve",
        &splinepy::py::ApproximateCurve,
        py::arg("points"),
        py::arg("degree"),
        py::arg("n_control_points"),
        py::arg("centripetal"),
        py::arg("knot_vector"));
  m.def("interpolate_surface",
        &splinepy::py::InterpolateSurface,
        py::arg("points"),
        py::arg("size_u"),
        py::arg("size_v"),
        py::arg("degree_u"),
        py::arg("degree_v"),
        py::arg("centripetal"));
  m.def("approximate_surface",
        &splinepy::py::ApproximateSurface,
        py::arg("points"),
        py::arg("num_points_u"),
        py::arg("num_points_v"),
        py::arg("size_u"),
        py::arg("size_v"),
        py::arg("degree_u"),
        py::arg("degree_v"),
        py::arg("centripetal"));
}

} // namespace splinepy::py::init
