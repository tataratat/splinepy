#include <splinepy/py/py_nurbs.hpp>

void init_nurbs6(py::module_& m) {
  // 6P (Pamametric dimension)
  add_nurbs_pyclass<6, 1>(m, "NURBS6P1D");
  add_nurbs_pyclass<6, 2>(m, "NURBS6P2D");
  add_nurbs_pyclass<6, 3>(m, "NURBS6P3D");
  add_nurbs_pyclass<6, 4>(m, "NURBS6P4D");
  add_nurbs_pyclass<6, 5>(m, "NURBS6P5D");
  add_nurbs_pyclass<6, 6>(m, "NURBS6P6D");
  add_nurbs_pyclass<6, 7>(m, "NURBS6P7D");
  add_nurbs_pyclass<6, 8>(m, "NURBS6P8D");
  add_nurbs_pyclass<6, 9>(m, "NURBS6P9D");
  add_nurbs_pyclass<6, 10>(m, "NURBS6P10D");
}
