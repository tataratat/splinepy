#include <splinepy/py/py_rational_bezier.hpp>

void init_rational_bezier1(py::module_ &m) {
    // 1P (Pamametric dimension)
    add_rational_bezier_pyclass<1, 1>(m, "RationalBezier1P1D");
    add_rational_bezier_pyclass<1, 2>(m, "RationalBezier1P2D");
    add_rational_bezier_pyclass<1, 3>(m, "RationalBezier1P3D");
    add_rational_bezier_pyclass<1, 4>(m, "RationalBezier1P4D");
    add_rational_bezier_pyclass<1, 5>(m, "RationalBezier1P5D");
    add_rational_bezier_pyclass<1, 6>(m, "RationalBezier1P6D");
    add_rational_bezier_pyclass<1, 7>(m, "RationalBezier1P7D");
    add_rational_bezier_pyclass<1, 8>(m, "RationalBezier1P8D");
    add_rational_bezier_pyclass<1, 9>(m, "RationalBezier1P9D");
    add_rational_bezier_pyclass<1, 10>(m,"RationalBezier1P10D");
}
