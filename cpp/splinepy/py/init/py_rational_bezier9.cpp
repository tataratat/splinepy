#include <splinepy/py/py_rational_bezier.hpp>

void init_rational_bezier9(py::module_& m) {
  // 9P (Pamametric dimension)
  add_rational_bezier_pyclass<9, 1>(m, "RationalBezier9P1D");
  add_rational_bezier_pyclass<9, 2>(m, "RationalBezier9P2D");
  add_rational_bezier_pyclass<9, 3>(m, "RationalBezier9P3D");
  add_rational_bezier_pyclass<9, 4>(m, "RationalBezier9P4D");
  add_rational_bezier_pyclass<9, 5>(m, "RationalBezier9P5D");
  add_rational_bezier_pyclass<9, 6>(m, "RationalBezier9P6D");
  add_rational_bezier_pyclass<9, 7>(m, "RationalBezier9P7D");
  add_rational_bezier_pyclass<9, 8>(m, "RationalBezier9P8D");
  add_rational_bezier_pyclass<9, 9>(m, "RationalBezier9P9D");
  add_rational_bezier_pyclass<9, 10>(m, "RationalBezier9P10D");
}
