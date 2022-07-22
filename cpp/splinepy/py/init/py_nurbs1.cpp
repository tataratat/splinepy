#include <splinepy/py/py_nurbs.hpp>

void init_nurbs1(py::module_ &m) {
    // 1P (Pamametric dimension)
    add_nurbs_pyclass<1, 1>(m, "NURBS1P1D");
    add_nurbs_pyclass<1, 2>(m, "NURBS1P2D");
    add_nurbs_pyclass<1, 3>(m, "NURBS1P3D");
    add_nurbs_pyclass<1, 4>(m, "NURBS1P4D");
    add_nurbs_pyclass<1, 5>(m, "NURBS1P5D");
    add_nurbs_pyclass<1, 6>(m, "NURBS1P6D");
    add_nurbs_pyclass<1, 7>(m, "NURBS1P7D");
    add_nurbs_pyclass<1, 8>(m, "NURBS1P8D");
    add_nurbs_pyclass<1, 9>(m, "NURBS1P9D");
    add_nurbs_pyclass<1, 10>(m, "NURBS1P10D");
}
