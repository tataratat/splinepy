#include <bspline.hpp>

void init_bspline1(py::module_ &m) {
    // 1P (Pamametric dimension)
    add_bspline_pyclass<1, 1>(m, "BSpline1P1D");
    add_bspline_pyclass<1, 2>(m, "BSpline1P2D");
    add_bspline_pyclass<1, 3>(m, "BSpline1P3D");
    add_bspline_pyclass<1, 4>(m, "BSpline1P4D");
    add_bspline_pyclass<1, 5>(m, "BSpline1P5D");
    add_bspline_pyclass<1, 6>(m, "BSpline1P6D");
    add_bspline_pyclass<1, 7>(m, "BSpline1P7D");
    add_bspline_pyclass<1, 8>(m, "BSpline1P8D");
    add_bspline_pyclass<1, 9>(m, "BSpline1P9D");
    add_bspline_pyclass<1, 10>(m, "BSpline1P10D");
}
