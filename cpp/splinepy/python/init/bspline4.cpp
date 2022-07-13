#include <bspline.hpp>

void init_bspline4(py::module_ &m) {
    // 4P (Pamametric dimension)
    add_bspline_pyclass<4, 1>(m, "BSpline4P1D");
    add_bspline_pyclass<4, 2>(m, "BSpline4P2D");
    add_bspline_pyclass<4, 3>(m, "BSpline4P3D");
    add_bspline_pyclass<4, 4>(m, "BSpline4P4D");
    add_bspline_pyclass<4, 5>(m, "BSpline4P5D");
    add_bspline_pyclass<4, 6>(m, "BSpline4P6D");
    add_bspline_pyclass<4, 7>(m, "BSpline4P7D");
    add_bspline_pyclass<4, 8>(m, "BSpline4P8D");
    add_bspline_pyclass<4, 9>(m, "BSpline4P9D");
    add_bspline_pyclass<4, 10>(m, "BSpline4P10D");
}
