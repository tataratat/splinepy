#include <splinepy/py/py_nurbs.hpp>

void init_nurbs7(py::module_ &m) {
    // 7P (Pamametric dimension)
    add_nurbs_pyclass<7, 1>(m, "NURBS7P1D");
    add_nurbs_pyclass<7, 2>(m, "NURBS7P2D");
    add_nurbs_pyclass<7, 3>(m, "NURBS7P3D");
    add_nurbs_pyclass<7, 4>(m, "NURBS7P4D");
    add_nurbs_pyclass<7, 5>(m, "NURBS7P5D");
    add_nurbs_pyclass<7, 6>(m, "NURBS7P6D");
    add_nurbs_pyclass<7, 7>(m, "NURBS7P7D");
    add_nurbs_pyclass<7, 8>(m, "NURBS7P8D");
    add_nurbs_pyclass<7, 9>(m, "NURBS7P9D");
    add_nurbs_pyclass<7, 10>(m, "NURBS7P10D");
}
