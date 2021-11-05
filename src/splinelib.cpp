#include <spline_reader.hpp>
#include <bspline.hpp>
#include <nurbs.hpp>

PYBIND11_MODULE(_splinelibpy, m) {

    /* BSPLINE */

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

    // 2P (Pamametric dimension)
    add_bspline_pyclass<2, 1>(m, "BSpline2P1D");
    add_bspline_pyclass<2, 2>(m, "BSpline2P2D");
    add_bspline_pyclass<2, 3>(m, "BSpline2P3D");
    add_bspline_pyclass<2, 4>(m, "BSpline2P4D");
    add_bspline_pyclass<2, 5>(m, "BSpline2P5D");
    add_bspline_pyclass<2, 6>(m, "BSpline2P6D");
    add_bspline_pyclass<2, 7>(m, "BSpline2P7D");
    add_bspline_pyclass<2, 8>(m, "BSpline2P8D");
    add_bspline_pyclass<2, 9>(m, "BSpline2P9D");
    add_bspline_pyclass<2, 10>(m, "BSpline2P10D");

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

    // 6P (Pamametric dimension)
    add_bspline_pyclass<6, 1>(m, "BSpline6P1D");
    add_bspline_pyclass<6, 2>(m, "BSpline6P2D");
    add_bspline_pyclass<6, 3>(m, "BSpline6P3D");
    add_bspline_pyclass<6, 4>(m, "BSpline6P4D");
    add_bspline_pyclass<6, 5>(m, "BSpline6P5D");
    add_bspline_pyclass<6, 6>(m, "BSpline6P6D");
    add_bspline_pyclass<6, 7>(m, "BSpline6P7D");
    add_bspline_pyclass<6, 8>(m, "BSpline6P8D");
    add_bspline_pyclass<6, 9>(m, "BSpline6P9D");
    add_bspline_pyclass<6, 10>(m, "BSpline6P10D");

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

    // 8P (Pamametric dimension)
    add_bspline_pyclass<8, 1>(m, "BSpline8P1D");
    add_bspline_pyclass<8, 2>(m, "BSpline8P2D");
    add_bspline_pyclass<8, 3>(m, "BSpline8P3D");
    add_bspline_pyclass<8, 4>(m, "BSpline8P4D");
    add_bspline_pyclass<8, 5>(m, "BSpline8P5D");
    add_bspline_pyclass<8, 6>(m, "BSpline8P6D");
    add_bspline_pyclass<8, 7>(m, "BSpline8P7D");
    add_bspline_pyclass<8, 8>(m, "BSpline8P8D");
    add_bspline_pyclass<8, 9>(m, "BSpline8P9D");
    add_bspline_pyclass<8, 10>(m, "BSpline8P10D");

    // 9P (Pamametric dimension)
    add_bspline_pyclass<9, 1>(m, "BSpline9P1D");
    add_bspline_pyclass<9, 2>(m, "BSpline9P2D");
    add_bspline_pyclass<9, 3>(m, "BSpline9P3D");
    add_bspline_pyclass<9, 4>(m, "BSpline9P4D");
    add_bspline_pyclass<9, 5>(m, "BSpline9P5D");
    add_bspline_pyclass<9, 6>(m, "BSpline9P6D");
    add_bspline_pyclass<9, 7>(m, "BSpline9P7D");
    add_bspline_pyclass<9, 8>(m, "BSpline9P8D");
    add_bspline_pyclass<9, 9>(m, "BSpline9P9D");
    add_bspline_pyclass<9, 10>(m, "BSpline9P10D");

    // 10P (Pamametric dimension)
    add_bspline_pyclass<10, 1>(m, "BSpline10P1D");
    add_bspline_pyclass<10, 2>(m, "BSpline10P2D");
    add_bspline_pyclass<10, 3>(m, "BSpline10P3D");
    add_bspline_pyclass<10, 4>(m, "BSpline10P4D");
    add_bspline_pyclass<10, 5>(m, "BSpline10P5D");
    add_bspline_pyclass<10, 6>(m, "BSpline10P6D");
    add_bspline_pyclass<10, 7>(m, "BSpline10P7D");
    add_bspline_pyclass<10, 8>(m, "BSpline10P8D");
    add_bspline_pyclass<10, 9>(m, "BSpline10P9D");
    add_bspline_pyclass<10, 10>(m, "BSpline10P10D");

    /* BSPLINE END */

    /* NURBS */

    // 1P (Pamametric dimension)
    add_nurbs_pyclass<1, 1>(m, "NURBS1P1D");
    add_nurbs_pyclass<1, 2>(m, "NURBS1P2D");
    add_nurbs_pyclass<1, 3>(m, "NURBS1P3D");
    add_nurbs_pyclass<1, 4>(m, "NURBS1P4D");
    add_nurbs_pyclass<1, 5>(m, "NURBS1P5D");
    add_nurbs_pyclass<1, 6>(m, "NURBS1P6D");
    add_nurbs_pyclass<1, 7>(m, "NURBS1P7D");
    add_nurbs_pyclass<1, 8>(m, "NURBS1P8D");
    add_nurbs_pyclass<1, 9>(m, "NURBS1P9D");
    add_nurbs_pyclass<1, 10>(m, "NURBS1P10D");

    // 2P (Pamametric dimension)
    add_nurbs_pyclass<2, 1>(m, "NURBS2P1D");
    add_nurbs_pyclass<2, 2>(m, "NURBS2P2D");
    add_nurbs_pyclass<2, 3>(m, "NURBS2P3D");
    add_nurbs_pyclass<2, 4>(m, "NURBS2P4D");
    add_nurbs_pyclass<2, 5>(m, "NURBS2P5D");
    add_nurbs_pyclass<2, 6>(m, "NURBS2P6D");
    add_nurbs_pyclass<2, 7>(m, "NURBS2P7D");
    add_nurbs_pyclass<2, 8>(m, "NURBS2P8D");
    add_nurbs_pyclass<2, 9>(m, "NURBS2P9D");
    add_nurbs_pyclass<2, 10>(m, "NURBS2P10D");

    // 3P (Pamametric dimension)
    add_nurbs_pyclass<3, 1>(m, "NURBS3P1D");
    add_nurbs_pyclass<3, 2>(m, "NURBS3P2D");
    add_nurbs_pyclass<3, 3>(m, "NURBS3P3D");
    add_nurbs_pyclass<3, 4>(m, "NURBS3P4D");
    add_nurbs_pyclass<3, 5>(m, "NURBS3P5D");
    add_nurbs_pyclass<3, 6>(m, "NURBS3P6D");
    add_nurbs_pyclass<3, 7>(m, "NURBS3P7D");
    add_nurbs_pyclass<3, 8>(m, "NURBS3P8D");
    add_nurbs_pyclass<3, 9>(m, "NURBS3P9D");
    add_nurbs_pyclass<3, 10>(m, "NURBS3P10D");

    // 4P (Pamametric dimension)
    add_nurbs_pyclass<4, 1>(m, "NURBS4P1D");
    add_nurbs_pyclass<4, 2>(m, "NURBS4P2D");
    add_nurbs_pyclass<4, 3>(m, "NURBS4P3D");
    add_nurbs_pyclass<4, 4>(m, "NURBS4P4D");
    add_nurbs_pyclass<4, 5>(m, "NURBS4P5D");
    add_nurbs_pyclass<4, 6>(m, "NURBS4P6D");
    add_nurbs_pyclass<4, 7>(m, "NURBS4P7D");
    add_nurbs_pyclass<4, 8>(m, "NURBS4P8D");
    add_nurbs_pyclass<4, 9>(m, "NURBS4P9D");
    add_nurbs_pyclass<4, 10>(m, "NURBS4P10D");

    // 5P (Pamametric dimension)
    add_nurbs_pyclass<5, 1>(m, "NURBS5P1D");
    add_nurbs_pyclass<5, 2>(m, "NURBS5P2D");
    add_nurbs_pyclass<5, 3>(m, "NURBS5P3D");
    add_nurbs_pyclass<5, 4>(m, "NURBS5P4D");
    add_nurbs_pyclass<5, 5>(m, "NURBS5P5D");
    add_nurbs_pyclass<5, 6>(m, "NURBS5P6D");
    add_nurbs_pyclass<5, 7>(m, "NURBS5P7D");
    add_nurbs_pyclass<5, 8>(m, "NURBS5P8D");
    add_nurbs_pyclass<5, 9>(m, "NURBS5P9D");
    add_nurbs_pyclass<5, 10>(m, "NURBS5P10D");

    // 6P (Pamametric dimension)
    add_nurbs_pyclass<6, 1>(m, "NURBS6P1D");
    add_nurbs_pyclass<6, 2>(m, "NURBS6P2D");
    add_nurbs_pyclass<6, 3>(m, "NURBS6P3D");
    add_nurbs_pyclass<6, 4>(m, "NURBS6P4D");
    add_nurbs_pyclass<6, 5>(m, "NURBS6P5D");
    add_nurbs_pyclass<6, 6>(m, "NURBS6P6D");
    add_nurbs_pyclass<6, 7>(m, "NURBS6P7D");
    add_nurbs_pyclass<6, 8>(m, "NURBS6P8D");
    add_nurbs_pyclass<6, 9>(m, "NURBS6P9D");
    add_nurbs_pyclass<6, 10>(m, "NURBS6P10D");

    // 7P (Pamametric dimension)
    add_nurbs_pyclass<7, 1>(m, "NURBS7P1D");
    add_nurbs_pyclass<7, 2>(m, "NURBS7P2D");
    add_nurbs_pyclass<7, 3>(m, "NURBS7P3D");
    add_nurbs_pyclass<7, 4>(m, "NURBS7P4D");
    add_nurbs_pyclass<7, 5>(m, "NURBS7P5D");
    add_nurbs_pyclass<7, 6>(m, "NURBS7P6D");
    add_nurbs_pyclass<7, 7>(m, "NURBS7P7D");
    add_nurbs_pyclass<7, 8>(m, "NURBS7P8D");
    add_nurbs_pyclass<7, 9>(m, "NURBS7P9D");
    add_nurbs_pyclass<7, 10>(m, "NURBS7P10D");

    // 8P (Pamametric dimension)
    add_nurbs_pyclass<8, 1>(m, "NURBS8P1D");
    add_nurbs_pyclass<8, 2>(m, "NURBS8P2D");
    add_nurbs_pyclass<8, 3>(m, "NURBS8P3D");
    add_nurbs_pyclass<8, 4>(m, "NURBS8P4D");
    add_nurbs_pyclass<8, 5>(m, "NURBS8P5D");
    add_nurbs_pyclass<8, 6>(m, "NURBS8P6D");
    add_nurbs_pyclass<8, 7>(m, "NURBS8P7D");
    add_nurbs_pyclass<8, 8>(m, "NURBS8P8D");
    add_nurbs_pyclass<8, 9>(m, "NURBS8P9D");
    add_nurbs_pyclass<8, 10>(m, "NURBS8P10D");

    // 9P (Pamametric dimension)
    add_nurbs_pyclass<9, 1>(m, "NURBS9P1D");
    add_nurbs_pyclass<9, 2>(m, "NURBS9P2D");
    add_nurbs_pyclass<9, 3>(m, "NURBS9P3D");
    add_nurbs_pyclass<9, 4>(m, "NURBS9P4D");
    add_nurbs_pyclass<9, 5>(m, "NURBS9P5D");
    add_nurbs_pyclass<9, 6>(m, "NURBS9P6D");
    add_nurbs_pyclass<9, 7>(m, "NURBS9P7D");
    add_nurbs_pyclass<9, 8>(m, "NURBS9P8D");
    add_nurbs_pyclass<9, 9>(m, "NURBS9P9D");
    add_nurbs_pyclass<9, 10>(m, "NURBS9P10D");

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

    /* NURBS END */

    /* SPLINE READER */
    py::class_<SplineReader>(m, "Reader")
        .def(py::init<>())
        .def("read_iges", &SplineReader::read_iges, py::arg("fname"))
        .def("read_xml", &SplineReader::read_xml, py::arg("fname"))
        .def("read_irit", &SplineReader::read_irit, py::arg("fname"))
    ;
    /* SPLINE READER END */

}
