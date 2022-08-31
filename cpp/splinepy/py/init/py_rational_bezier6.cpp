#include <splinepy/py/py_rational_bezier.hpp>

void init_rational_bezier6(py::module_& m) {
  // 6P (Pamametric dimension)
  add_rational_bezier_pyclass<6, 1>(m, "RationalBezier6P1D");
  add_rational_bezier_pyclass<6, 2>(m, "RationalBezier6P2D");
  add_rational_bezier_pyclass<6, 3>(m, "RationalBezier6P3D");
  add_rational_bezier_pyclass<6, 4>(m, "RationalBezier6P4D");
  add_rational_bezier_pyclass<6, 5>(m, "RationalBezier6P5D");
  add_rational_bezier_pyclass<6, 6>(m, "RationalBezier6P6D");
  add_rational_bezier_pyclass<6, 7>(m, "RationalBezier6P7D");
  add_rational_bezier_pyclass<6, 8>(m, "RationalBezier6P8D");
  add_rational_bezier_pyclass<6, 9>(m, "RationalBezier6P9D");
  add_rational_bezier_pyclass<6, 10>(m, "RationalBezier6P10D");
}
