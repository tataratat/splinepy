#include <splinepy/py/py_bspline.hpp>

void init_bspline2(py::module_& m) {
  // 2P (Pamametric dimension)
  add_bspline_pyclass<2, 1>(m, "BSpline2P1D");
  add_bspline_pyclass<2, 2>(m, "BSpline2P2D");
  add_bspline_pyclass<2, 3>(m, "BSpline2P3D");
  add_bspline_pyclass<2, 4>(m, "BSpline2P4D");
  add_bspline_pyclass<2, 5>(m, "BSpline2P5D");
  add_bspline_pyclass<2, 6>(m, "BSpline2P6D");
  add_bspline_pyclass<2, 7>(m, "BSpline2P7D");
  add_bspline_pyclass<2, 8>(m, "BSpline2P8D");
  add_bspline_pyclass<2, 9>(m, "BSpline2P9D");
  add_bspline_pyclass<2, 10>(m, "BSpline2P10D");
}
