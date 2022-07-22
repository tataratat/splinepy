#include <splinepy/py/py_bezier.hpp>

void init_bezier1(py::module_ &m) {
    // 1P (Pamametric dimension)
    add_bezier_pyclass<1, 1>(m, "Bezier1P1D");
    add_bezier_pyclass<1, 2>(m, "Bezier1P2D");
    add_bezier_pyclass<1, 3>(m, "Bezier1P3D");
    add_bezier_pyclass<1, 4>(m, "Bezier1P4D");
    add_bezier_pyclass<1, 5>(m, "Bezier1P5D");
    add_bezier_pyclass<1, 6>(m, "Bezier1P6D");
    add_bezier_pyclass<1, 7>(m, "Bezier1P7D");
    add_bezier_pyclass<1, 8>(m, "Bezier1P8D");
    add_bezier_pyclass<1, 9>(m, "Bezier1P9D");
    add_bezier_pyclass<1, 10>(m, "Bezier1P10D");
}
