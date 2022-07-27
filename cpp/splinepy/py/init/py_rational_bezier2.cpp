#include <splinepy/py/py_rational_bezier.hpp>

void init_rational_bezier2(py::module_ &m) {
    // 2P (Pamametric dimension)
    add_rational_bezier_pyclass<2, 1>(m, "RationalBezier2P1D");
    add_rational_bezier_pyclass<2, 2>(m, "RationalBezier2P2D");
    add_rational_bezier_pyclass<2, 3>(m, "RationalBezier2P3D");
    add_rational_bezier_pyclass<2, 4>(m, "RationalBezier2P4D");
    add_rational_bezier_pyclass<2, 5>(m, "RationalBezier2P5D");
    add_rational_bezier_pyclass<2, 6>(m, "RationalBezier2P6D");
    add_rational_bezier_pyclass<2, 7>(m, "RationalBezier2P7D");
    add_rational_bezier_pyclass<2, 8>(m, "RationalBezier2P8D");
    add_rational_bezier_pyclass<2, 9>(m, "RationalBezier2P9D");
    add_rational_bezier_pyclass<2, 10>(m,"RationalBezier2P10D");
}
