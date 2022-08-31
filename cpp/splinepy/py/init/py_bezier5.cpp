#include <splinepy/py/py_bezier.hpp>
#include <splinepy/py/py_rational_bezier.hpp>

void init_bezier5(py::module_& m) {
  // 5P (Pamametric dimension)
  add_bezier_pyclass<5, 1>(m, "Bezier5P1D");
  add_bezier_pyclass<5, 2>(m, "Bezier5P2D");
  add_bezier_pyclass<5, 3>(m, "Bezier5P3D");
  add_bezier_pyclass<5, 4>(m, "Bezier5P4D");
  add_bezier_pyclass<5, 5>(m, "Bezier5P5D");
  add_bezier_pyclass<5, 6>(m, "Bezier5P6D");
  add_bezier_pyclass<5, 7>(m, "Bezier5P7D");
  add_bezier_pyclass<5, 8>(m, "Bezier5P8D");
  add_bezier_pyclass<5, 9>(m, "Bezier5P9D");
  add_bezier_pyclass<5, 10>(m, "Bezier5P10D");
}
