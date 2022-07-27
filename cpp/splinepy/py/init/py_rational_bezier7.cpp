#include <splinepy/py/py_rational_bezier.hpp>

void init_rational_bezier7(py::module_ &m) {
    // 7P (Pamametric dimension)
    add_rational_bezier_pyclass<7, 1>(m, "RationalBezier7P1D");
    add_rational_bezier_pyclass<7, 2>(m, "RationalBezier7P2D");
    add_rational_bezier_pyclass<7, 3>(m, "RationalBezier7P3D");
    add_rational_bezier_pyclass<7, 4>(m, "RationalBezier7P4D");
    add_rational_bezier_pyclass<7, 5>(m, "RationalBezier7P5D");
    add_rational_bezier_pyclass<7, 6>(m, "RationalBezier7P6D");
    add_rational_bezier_pyclass<7, 7>(m, "RationalBezier7P7D");
    add_rational_bezier_pyclass<7, 8>(m, "RationalBezier7P8D");
    add_rational_bezier_pyclass<7, 9>(m, "RationalBezier7P9D");
    add_rational_bezier_pyclass<7, 10>(m,"RationalBezier7P10D");
}
