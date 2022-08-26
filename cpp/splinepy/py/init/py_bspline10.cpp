#include <splinepy/py/py_bspline.hpp>

void init_bspline10(py::module_& m) {
  // 10P (Pamametric dimension)
  add_bspline_pyclass<10, 1>(m, "BSpline10P1D");
  add_bspline_pyclass<10, 2>(m, "BSpline10P2D");
  add_bspline_pyclass<10, 3>(m, "BSpline10P3D");
  add_bspline_pyclass<10, 4>(m, "BSpline10P4D");
  add_bspline_pyclass<10, 5>(m, "BSpline10P5D");
  add_bspline_pyclass<10, 6>(m, "BSpline10P6D");
  add_bspline_pyclass<10, 7>(m, "BSpline10P7D");
  add_bspline_pyclass<10, 8>(m, "BSpline10P8D");
  add_bspline_pyclass<10, 9>(m, "BSpline10P9D");
  add_bspline_pyclass<10, 10>(m, "BSpline10P10D");
}
