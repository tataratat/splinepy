#include <splinepy/py/py_bezier.hpp>
#include <splinepy/py/py_rational_bezier.hpp>

void init_bezier2(py::module_ &m) {
    // 2P (Pamametric dimension)
    add_bezier_pyclass<2, 1>(m, "Bezier2P1D");
    add_bezier_pyclass<2, 2>(m, "Bezier2P2D");
    add_bezier_pyclass<2, 3>(m, "Bezier2P3D");
    add_bezier_pyclass<2, 4>(m, "Bezier2P4D");
    add_bezier_pyclass<2, 5>(m, "Bezier2P5D");
    add_bezier_pyclass<2, 6>(m, "Bezier2P6D");
    add_bezier_pyclass<2, 7>(m, "Bezier2P7D");
    add_bezier_pyclass<2, 8>(m, "Bezier2P8D");
    add_bezier_pyclass<2, 9>(m, "Bezier2P9D");
    add_bezier_pyclass<2, 10>(m, "Bezier2P10D");
}
