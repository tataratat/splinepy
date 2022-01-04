#include <nurbs.hpp>

void init_nurbs8(py::module_ &m) {
    // 8P (Pamametric dimension)
    add_nurbs_pyclass<8, 1>(m, "NURBS8P1D");
    add_nurbs_pyclass<8, 2>(m, "NURBS8P2D");
    add_nurbs_pyclass<8, 3>(m, "NURBS8P3D");
    add_nurbs_pyclass<8, 4>(m, "NURBS8P4D");
    add_nurbs_pyclass<8, 5>(m, "NURBS8P5D");
    add_nurbs_pyclass<8, 6>(m, "NURBS8P6D");
    add_nurbs_pyclass<8, 7>(m, "NURBS8P7D");
    add_nurbs_pyclass<8, 8>(m, "NURBS8P8D");
    add_nurbs_pyclass<8, 9>(m, "NURBS8P9D");
    add_nurbs_pyclass<8, 10>(m, "NURBS8P10D");
}
