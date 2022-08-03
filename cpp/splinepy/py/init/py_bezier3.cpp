#include <splinepy/py/py_bezier.hpp>
#include <splinepy/py/py_rational_bezier.hpp>

void init_bezier3(py::module_ &m) {
    // 3P (Pamametric dimension)
    add_bezier_pyclass<3, 1>(m, "Bezier3P1D");
    add_bezier_pyclass<3, 2>(m, "Bezier3P2D");
    add_bezier_pyclass<3, 3>(m, "Bezier3P3D");
    add_bezier_pyclass<3, 4>(m, "Bezier3P4D");
    add_bezier_pyclass<3, 5>(m, "Bezier3P5D");
    add_bezier_pyclass<3, 6>(m, "Bezier3P6D");
    add_bezier_pyclass<3, 7>(m, "Bezier3P7D");
    add_bezier_pyclass<3, 8>(m, "Bezier3P8D");
    add_bezier_pyclass<3, 9>(m, "Bezier3P9D");
    add_bezier_pyclass<3, 10>(m, "Bezier3P10D");
}
