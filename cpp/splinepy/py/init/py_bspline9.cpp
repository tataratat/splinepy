#include <splinepy/py/py_bspline.hpp>

void init_bspline9(py::module_& m) {
  // 9P (Pamametric dimension)
  add_bspline_pyclass<9, 1>(m, "BSpline9P1D");
  add_bspline_pyclass<9, 2>(m, "BSpline9P2D");
  add_bspline_pyclass<9, 3>(m, "BSpline9P3D");
  add_bspline_pyclass<9, 4>(m, "BSpline9P4D");
  add_bspline_pyclass<9, 5>(m, "BSpline9P5D");
  add_bspline_pyclass<9, 6>(m, "BSpline9P6D");
  add_bspline_pyclass<9, 7>(m, "BSpline9P7D");
  add_bspline_pyclass<9, 8>(m, "BSpline9P8D");
  add_bspline_pyclass<9, 9>(m, "BSpline9P9D");
  add_bspline_pyclass<9, 10>(m, "BSpline9P10D");
}
