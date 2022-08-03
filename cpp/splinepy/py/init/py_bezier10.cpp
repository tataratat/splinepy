#include <splinepy/py/py_bezier.hpp>
#include <splinepy/py/py_rational_bezier.hpp>

void init_bezier10(py::module_ &m) {
    // 10P (Pamametric dimension)
    add_bezier_pyclass<10, 1>(m, "Bezier10P1D");
    add_bezier_pyclass<10, 2>(m, "Bezier10P2D");
    add_bezier_pyclass<10, 3>(m, "Bezier10P3D");
    add_bezier_pyclass<10, 4>(m, "Bezier10P4D");
    add_bezier_pyclass<10, 5>(m, "Bezier10P5D");
    add_bezier_pyclass<10, 6>(m, "Bezier10P6D");
    add_bezier_pyclass<10, 7>(m, "Bezier10P7D");
    add_bezier_pyclass<10, 8>(m, "Bezier10P8D");
    add_bezier_pyclass<10, 9>(m, "Bezier10P9D");
    add_bezier_pyclass<10, 10>(m, "Bezier10P10D");
}
