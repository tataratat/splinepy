#include <splinepy/py/py_rational_bezier.hpp>

void init_rational_bezier8(py::module_& m) {
  // 8P (Pamametric dimension)
  add_rational_bezier_pyclass<8, 1>(m, "RationalBezier8P1D");
  add_rational_bezier_pyclass<8, 2>(m, "RationalBezier8P2D");
  add_rational_bezier_pyclass<8, 3>(m, "RationalBezier8P3D");
  add_rational_bezier_pyclass<8, 4>(m, "RationalBezier8P4D");
  add_rational_bezier_pyclass<8, 5>(m, "RationalBezier8P5D");
  add_rational_bezier_pyclass<8, 6>(m, "RationalBezier8P6D");
  add_rational_bezier_pyclass<8, 7>(m, "RationalBezier8P7D");
  add_rational_bezier_pyclass<8, 8>(m, "RationalBezier8P8D");
  add_rational_bezier_pyclass<8, 9>(m, "RationalBezier8P9D");
  add_rational_bezier_pyclass<8, 10>(m, "RationalBezier8P10D");
}
