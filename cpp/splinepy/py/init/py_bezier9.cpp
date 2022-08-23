#include <splinepy/py/py_bezier.hpp>
#include <splinepy/py/py_rational_bezier.hpp>

void init_bezier9(py::module_& m) {
  // 9P (Pamametric dimension)
  add_bezier_pyclass<9, 1>(m, "Bezier9P1D");
  add_bezier_pyclass<9, 2>(m, "Bezier9P2D");
  add_bezier_pyclass<9, 3>(m, "Bezier9P3D");
  add_bezier_pyclass<9, 4>(m, "Bezier9P4D");
  add_bezier_pyclass<9, 5>(m, "Bezier9P5D");
  add_bezier_pyclass<9, 6>(m, "Bezier9P6D");
  add_bezier_pyclass<9, 7>(m, "Bezier9P7D");
  add_bezier_pyclass<9, 8>(m, "Bezier9P8D");
  add_bezier_pyclass<9, 9>(m, "Bezier9P9D");
  add_bezier_pyclass<9, 10>(m, "Bezier9P10D");
}
