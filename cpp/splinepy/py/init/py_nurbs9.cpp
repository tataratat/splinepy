#include <splinepy/py/py_nurbs.hpp>

void init_nurbs9(py::module_& m) {
  // 9P (Pamametric dimension)
  add_nurbs_pyclass<9, 1>(m, "NURBS9P1D");
  add_nurbs_pyclass<9, 2>(m, "NURBS9P2D");
  add_nurbs_pyclass<9, 3>(m, "NURBS9P3D");
  add_nurbs_pyclass<9, 4>(m, "NURBS9P4D");
  add_nurbs_pyclass<9, 5>(m, "NURBS9P5D");
  add_nurbs_pyclass<9, 6>(m, "NURBS9P6D");
  add_nurbs_pyclass<9, 7>(m, "NURBS9P7D");
  add_nurbs_pyclass<9, 8>(m, "NURBS9P8D");
  add_nurbs_pyclass<9, 9>(m, "NURBS9P9D");
  add_nurbs_pyclass<9, 10>(m, "NURBS9P10D");
}
