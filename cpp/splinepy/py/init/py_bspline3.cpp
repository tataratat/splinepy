#include <splinepy/py/py_bspline.hpp>

void init_bspline3(py::module_ &m) {
    // 3P (Pamametric dimension)
    add_bspline_pyclass<3, 1>(m, "BSpline3P1D");
    add_bspline_pyclass<3, 2>(m, "BSpline3P2D");
    add_bspline_pyclass<3, 3>(m, "BSpline3P3D");
    add_bspline_pyclass<3, 4>(m, "BSpline3P4D");
    add_bspline_pyclass<3, 5>(m, "BSpline3P5D");
    add_bspline_pyclass<3, 6>(m, "BSpline3P6D");
    add_bspline_pyclass<3, 7>(m, "BSpline3P7D");
    add_bspline_pyclass<3, 8>(m, "BSpline3P8D");
    add_bspline_pyclass<3, 9>(m, "BSpline3P9D");
    add_bspline_pyclass<3, 10>(m, "BSpline3P10D");
}
