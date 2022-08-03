#include <splinepy/py/py_bezier.hpp>
#include <splinepy/py/py_rational_bezier.hpp>

void init_bezier6(py::module_ &m) {
    // 6P (Pamametric dimension)
    add_bezier_pyclass<6, 1>(m, "Bezier6P1D");
    add_bezier_pyclass<6, 2>(m, "Bezier6P2D");
    add_bezier_pyclass<6, 3>(m, "Bezier6P3D");
    add_bezier_pyclass<6, 4>(m, "Bezier6P4D");
    add_bezier_pyclass<6, 5>(m, "Bezier6P5D");
    add_bezier_pyclass<6, 6>(m, "Bezier6P6D");
    add_bezier_pyclass<6, 7>(m, "Bezier6P7D");
    add_bezier_pyclass<6, 8>(m, "Bezier6P8D");
    add_bezier_pyclass<6, 9>(m, "Bezier6P9D");
    add_bezier_pyclass<6, 10>(m, "Bezier6P10D");
}
