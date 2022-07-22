#include <splinepy/py/py_bezier.hpp>

void init_bezier7(py::module_ &m) {
    // 7P (Pamametric dimension)
    add_bezier_pyclass<7, 1>(m, "Bezier7P1D");
    add_bezier_pyclass<7, 2>(m, "Bezier7P2D");
    add_bezier_pyclass<7, 3>(m, "Bezier7P3D");
    add_bezier_pyclass<7, 4>(m, "Bezier7P4D");
    add_bezier_pyclass<7, 5>(m, "Bezier7P5D");
    add_bezier_pyclass<7, 6>(m, "Bezier7P6D");
    add_bezier_pyclass<7, 7>(m, "Bezier7P7D");
    add_bezier_pyclass<7, 8>(m, "Bezier7P8D");
    add_bezier_pyclass<7, 9>(m, "Bezier7P9D");
    add_bezier_pyclass<7, 10>(m, "Bezier7P10D");
}
