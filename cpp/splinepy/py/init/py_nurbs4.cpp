#include <splinepy/py/py_nurbs.hpp>

void init_nurbs4(py::module_ &m) {
    // 4P (Pamametric dimension)
    add_nurbs_pyclass<4, 1>(m, "NURBS4P1D");
    add_nurbs_pyclass<4, 2>(m, "NURBS4P2D");
    add_nurbs_pyclass<4, 3>(m, "NURBS4P3D");
    add_nurbs_pyclass<4, 4>(m, "NURBS4P4D");
    add_nurbs_pyclass<4, 5>(m, "NURBS4P5D");
    add_nurbs_pyclass<4, 6>(m, "NURBS4P6D");
    add_nurbs_pyclass<4, 7>(m, "NURBS4P7D");
    add_nurbs_pyclass<4, 8>(m, "NURBS4P8D");
    add_nurbs_pyclass<4, 9>(m, "NURBS4P9D");
    add_nurbs_pyclass<4, 10>(m, "NURBS4P10D");
}
