#include <nurbs.hpp>

void init_nurbs5(py::module_ &m) {
    // 5P (Pamametric dimension)
    add_nurbs_pyclass<5, 1>(m, "NURBS5P1D");
    add_nurbs_pyclass<5, 2>(m, "NURBS5P2D");
    add_nurbs_pyclass<5, 3>(m, "NURBS5P3D");
    add_nurbs_pyclass<5, 4>(m, "NURBS5P4D");
    add_nurbs_pyclass<5, 5>(m, "NURBS5P5D");
    add_nurbs_pyclass<5, 6>(m, "NURBS5P6D");
    add_nurbs_pyclass<5, 7>(m, "NURBS5P7D");
    add_nurbs_pyclass<5, 8>(m, "NURBS5P8D");
    add_nurbs_pyclass<5, 9>(m, "NURBS5P9D");
    add_nurbs_pyclass<5, 10>(m, "NURBS5P10D");
}
