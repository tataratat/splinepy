#include <splinepy/py/py_bezier.hpp>

void init_bezier4(py::module_ &m) {
    // 4P (Pamametric dimension)
    add_bezier_pyclass<4, 1>(m, "Bezier4P1D");
    add_bezier_pyclass<4, 2>(m, "Bezier4P2D");
    add_bezier_pyclass<4, 3>(m, "Bezier4P3D");
    add_bezier_pyclass<4, 4>(m, "Bezier4P4D");
    add_bezier_pyclass<4, 5>(m, "Bezier4P5D");
    add_bezier_pyclass<4, 6>(m, "Bezier4P6D");
    add_bezier_pyclass<4, 7>(m, "Bezier4P7D");
    add_bezier_pyclass<4, 8>(m, "Bezier4P8D");
    add_bezier_pyclass<4, 9>(m, "Bezier4P9D");
    add_bezier_pyclass<4, 10>(m, "Bezier4P10D");
}
