#include <splinepy/py/py_rational_bezier.hpp>

void init_rational_bezier3(py::module_ &m) {
    // 3P (Pamametric dimension)
    add_rational_bezier_pyclass<3, 1>(m, "RationalBezier3P1D");
    add_rational_bezier_pyclass<3, 2>(m, "RationalBezier3P2D");
    add_rational_bezier_pyclass<3, 3>(m, "RationalBezier3P3D");
    add_rational_bezier_pyclass<3, 4>(m, "RationalBezier3P4D");
    add_rational_bezier_pyclass<3, 5>(m, "RationalBezier3P5D");
    add_rational_bezier_pyclass<3, 6>(m, "RationalBezier3P6D");
    add_rational_bezier_pyclass<3, 7>(m, "RationalBezier3P7D");
    add_rational_bezier_pyclass<3, 8>(m, "RationalBezier3P8D");
    add_rational_bezier_pyclass<3, 9>(m, "RationalBezier3P9D");
    add_rational_bezier_pyclass<3, 10>(m,"RationalBezier3P10D");
}
