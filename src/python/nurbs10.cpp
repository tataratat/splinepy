#include <nurbs.hpp>

void init_nurbs10(py::module_ &m) {
    // 10P (Pamametric dimension)
    add_nurbs_pyclass<10, 1>(m, "NURBS10P1D");
    add_nurbs_pyclass<10, 2>(m, "NURBS10P2D");
    add_nurbs_pyclass<10, 3>(m, "NURBS10P3D");
    add_nurbs_pyclass<10, 4>(m, "NURBS10P4D");
    add_nurbs_pyclass<10, 5>(m, "NURBS10P5D");
    add_nurbs_pyclass<10, 6>(m, "NURBS10P6D");
    add_nurbs_pyclass<10, 7>(m, "NURBS10P7D");
    add_nurbs_pyclass<10, 8>(m, "NURBS10P8D");
    add_nurbs_pyclass<10, 9>(m, "NURBS10P9D");
    add_nurbs_pyclass<10, 10>(m, "NURBS10P10D");
}
