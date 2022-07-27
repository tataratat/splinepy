#include <splinepy/py/py_rational_bezier.hpp>

void init_rational_bezier5(py::module_ &m) {
    // 5P (Pamametric dimension)
    add_rational_bezier_pyclass<5, 1>(m, "RationalBezier5P1D");
    add_rational_bezier_pyclass<5, 2>(m, "RationalBezier5P2D");
    add_rational_bezier_pyclass<5, 3>(m, "RationalBezier5P3D");
    add_rational_bezier_pyclass<5, 4>(m, "RationalBezier5P4D");
    add_rational_bezier_pyclass<5, 5>(m, "RationalBezier5P5D");
    add_rational_bezier_pyclass<5, 6>(m, "RationalBezier5P6D");
    add_rational_bezier_pyclass<5, 7>(m, "RationalBezier5P7D");
    add_rational_bezier_pyclass<5, 8>(m, "RationalBezier5P8D");
    add_rational_bezier_pyclass<5, 9>(m, "RationalBezier5P9D");
    add_rational_bezier_pyclass<5, 10>(m,"RationalBezier5P10D");
}
