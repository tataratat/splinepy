#include <splinepy/py/py_bspline.hpp>

void init_bspline8(py::module_& m) {
  // 8P (Pamametric dimension)
  add_bspline_pyclass<8, 1>(m, "BSpline8P1D");
  add_bspline_pyclass<8, 2>(m, "BSpline8P2D");
  add_bspline_pyclass<8, 3>(m, "BSpline8P3D");
  add_bspline_pyclass<8, 4>(m, "BSpline8P4D");
  add_bspline_pyclass<8, 5>(m, "BSpline8P5D");
  add_bspline_pyclass<8, 6>(m, "BSpline8P6D");
  add_bspline_pyclass<8, 7>(m, "BSpline8P7D");
  add_bspline_pyclass<8, 8>(m, "BSpline8P8D");
  add_bspline_pyclass<8, 9>(m, "BSpline8P9D");
  add_bspline_pyclass<8, 10>(m, "BSpline8P10D");
}
