#include <splinepy/py/py_rational_bezier.hpp>

void init_rational_bezier10(py::module_ &m) {
    // 10P (Pamametric dimension)
    add_rational_bezier_pyclass<10, 1>(m, "RationalBezier10P1D");
    add_rational_bezier_pyclass<10, 2>(m, "RationalBezier10P2D");
    add_rational_bezier_pyclass<10, 3>(m, "RationalBezier10P3D");
    add_rational_bezier_pyclass<10, 4>(m, "RationalBezier10P4D");
    add_rational_bezier_pyclass<10, 5>(m, "RationalBezier10P5D");
    add_rational_bezier_pyclass<10, 6>(m, "RationalBezier10P6D");
    add_rational_bezier_pyclass<10, 7>(m, "RationalBezier10P7D");
    add_rational_bezier_pyclass<10, 8>(m, "RationalBezier10P8D");
    add_rational_bezier_pyclass<10, 9>(m, "RationalBezier10P9D");
    add_rational_bezier_pyclass<10, 10>(m,"RationalBezier10P10D");
}
