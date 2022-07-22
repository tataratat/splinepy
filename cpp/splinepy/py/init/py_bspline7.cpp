#include <splinepy/py/py_bspline.hpp>

void init_bspline7(py::module_ &m) {
    // 7P (Pamametric dimension)
    add_bspline_pyclass<7, 1>(m, "BSpline7P1D");
    add_bspline_pyclass<7, 2>(m, "BSpline7P2D");
    add_bspline_pyclass<7, 3>(m, "BSpline7P3D");
    add_bspline_pyclass<7, 4>(m, "BSpline7P4D");
    add_bspline_pyclass<7, 5>(m, "BSpline7P5D");
    add_bspline_pyclass<7, 6>(m, "BSpline7P6D");
    add_bspline_pyclass<7, 7>(m, "BSpline7P7D");
    add_bspline_pyclass<7, 8>(m, "BSpline7P8D");
    add_bspline_pyclass<7, 9>(m, "BSpline7P9D");
    add_bspline_pyclass<7, 10>(m, "BSpline7P10D");
}
