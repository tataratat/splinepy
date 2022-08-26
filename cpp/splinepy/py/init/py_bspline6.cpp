#include <splinepy/py/py_bspline.hpp>

void init_bspline6(py::module_& m) {
  // 6P (Pamametric dimension)
  add_bspline_pyclass<6, 1>(m, "BSpline6P1D");
  add_bspline_pyclass<6, 2>(m, "BSpline6P2D");
  add_bspline_pyclass<6, 3>(m, "BSpline6P3D");
  add_bspline_pyclass<6, 4>(m, "BSpline6P4D");
  add_bspline_pyclass<6, 5>(m, "BSpline6P5D");
  add_bspline_pyclass<6, 6>(m, "BSpline6P6D");
  add_bspline_pyclass<6, 7>(m, "BSpline6P7D");
  add_bspline_pyclass<6, 8>(m, "BSpline6P8D");
  add_bspline_pyclass<6, 9>(m, "BSpline6P9D");
  add_bspline_pyclass<6, 10>(m, "BSpline6P10D");
}
