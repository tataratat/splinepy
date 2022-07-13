#include <nurbs.hpp>

void init_nurbs3(py::module_ &m) {
    // 3P (Pamametric dimension)
    add_nurbs_pyclass<3, 1>(m, "NURBS3P1D");
    add_nurbs_pyclass<3, 2>(m, "NURBS3P2D");
    add_nurbs_pyclass<3, 3>(m, "NURBS3P3D");
    add_nurbs_pyclass<3, 4>(m, "NURBS3P4D");
    add_nurbs_pyclass<3, 5>(m, "NURBS3P5D");
    add_nurbs_pyclass<3, 6>(m, "NURBS3P6D");
    add_nurbs_pyclass<3, 7>(m, "NURBS3P7D");
    add_nurbs_pyclass<3, 8>(m, "NURBS3P8D");
    add_nurbs_pyclass<3, 9>(m, "NURBS3P9D");
    add_nurbs_pyclass<3, 10>(m, "NURBS3P10D");
}
