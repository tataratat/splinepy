#include <splinepy/py/py_rational_bezier.hpp>

void init_rational_bezier4(py::module_ &m) {
    // 4P (Pamametric dimension)
    add_rational_bezier_pyclass<4, 1>(m, "RationalBezier4P1D");
    add_rational_bezier_pyclass<4, 2>(m, "RationalBezier4P2D");
    add_rational_bezier_pyclass<4, 3>(m, "RationalBezier4P3D");
    add_rational_bezier_pyclass<4, 4>(m, "RationalBezier4P4D");
    add_rational_bezier_pyclass<4, 5>(m, "RationalBezier4P5D");
    add_rational_bezier_pyclass<4, 6>(m, "RationalBezier4P6D");
    add_rational_bezier_pyclass<4, 7>(m, "RationalBezier4P7D");
    add_rational_bezier_pyclass<4, 8>(m, "RationalBezier4P8D");
    add_rational_bezier_pyclass<4, 9>(m, "RationalBezier4P9D");
    add_rational_bezier_pyclass<4, 10>(m,"RationalBezier4P10D");
}
