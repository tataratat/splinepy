#include <splinepy/py/py_nurbs.hpp>

void init_nurbs2(py::module_& m) {
  // 2P (Pamametric dimension)
  add_nurbs_pyclass<2, 1>(m, "NURBS2P1D");
  add_nurbs_pyclass<2, 2>(m, "NURBS2P2D");
  add_nurbs_pyclass<2, 3>(m, "NURBS2P3D");
  add_nurbs_pyclass<2, 4>(m, "NURBS2P4D");
  add_nurbs_pyclass<2, 5>(m, "NURBS2P5D");
  add_nurbs_pyclass<2, 6>(m, "NURBS2P6D");
  add_nurbs_pyclass<2, 7>(m, "NURBS2P7D");
  add_nurbs_pyclass<2, 8>(m, "NURBS2P8D");
  add_nurbs_pyclass<2, 9>(m, "NURBS2P9D");
  add_nurbs_pyclass<2, 10>(m, "NURBS2P10D");
}
