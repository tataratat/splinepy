#include <splinepy/py/py_bezier.hpp>

void init_bezier8(py::module_ &m) {
    // 8P (Pamametric dimension)
    add_bezier_pyclass<8, 1>(m, "Bezier8P1D");
    add_bezier_pyclass<8, 2>(m, "Bezier8P2D");
    add_bezier_pyclass<8, 3>(m, "Bezier8P3D");
    add_bezier_pyclass<8, 4>(m, "Bezier8P4D");
    add_bezier_pyclass<8, 5>(m, "Bezier8P5D");
    add_bezier_pyclass<8, 6>(m, "Bezier8P6D");
    add_bezier_pyclass<8, 7>(m, "Bezier8P7D");
    add_bezier_pyclass<8, 8>(m, "Bezier8P8D");
    add_bezier_pyclass<8, 9>(m, "Bezier8P9D");
    add_bezier_pyclass<8, 10>(m, "Bezier8P10D");
}
