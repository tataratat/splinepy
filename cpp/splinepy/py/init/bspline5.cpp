#include <bspline.hpp>

void init_bspline5(py::module_ &m) {
    // 5P (Pamametric dimension)
    add_bspline_pyclass<5, 1>(m, "BSpline5P1D");
    add_bspline_pyclass<5, 2>(m, "BSpline5P2D");
    add_bspline_pyclass<5, 3>(m, "BSpline5P3D");
    add_bspline_pyclass<5, 4>(m, "BSpline5P4D");
    add_bspline_pyclass<5, 5>(m, "BSpline5P5D");
    add_bspline_pyclass<5, 6>(m, "BSpline5P6D");
    add_bspline_pyclass<5, 7>(m, "BSpline5P7D");
    add_bspline_pyclass<5, 8>(m, "BSpline5P8D");
    add_bspline_pyclass<5, 9>(m, "BSpline5P9D");
    add_bspline_pyclass<5, 10>(m, "BSpline5P10D");
}
