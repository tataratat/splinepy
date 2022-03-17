#include <bspline.hpp>
#include <nurbs.hpp>

// only compile a minimal set of splines
// 4P X 4D
void init_minimal(py::module_ &m) {
  // BSpline
  add_bspline_pyclass<1, 1>(m, "BSpline1P1D");
  add_bspline_pyclass<1, 2>(m, "BSpline1P2D");
  add_bspline_pyclass<1, 3>(m, "BSpline1P3D");
  add_bspline_pyclass<1, 4>(m, "BSpline1P4D");

  add_bspline_pyclass<2, 1>(m, "BSpline2P1D");
  add_bspline_pyclass<2, 2>(m, "BSpline2P2D");
  add_bspline_pyclass<2, 3>(m, "BSpline2P3D");
  add_bspline_pyclass<2, 4>(m, "BSpline2P4D");

  add_bspline_pyclass<3, 1>(m, "BSpline3P1D");
  add_bspline_pyclass<3, 2>(m, "BSpline3P2D");
  add_bspline_pyclass<3, 3>(m, "BSpline3P3D");
  add_bspline_pyclass<3, 4>(m, "BSpline3P4D");

  add_bspline_pyclass<4, 1>(m, "BSpline4P1D");
  add_bspline_pyclass<4, 2>(m, "BSpline4P2D");
  add_bspline_pyclass<4, 3>(m, "BSpline4P3D");
  add_bspline_pyclass<4, 4>(m, "BSpline4P4D");

  // NURBS
  add_nurbs_pyclass<1, 1>(m, "NURBS1P1D");
  add_nurbs_pyclass<1, 2>(m, "NURBS1P2D");
  add_nurbs_pyclass<1, 3>(m, "NURBS1P3D");
  add_nurbs_pyclass<1, 4>(m, "NURBS1P4D");

  add_nurbs_pyclass<2, 1>(m, "NURBS2P1D");
  add_nurbs_pyclass<2, 2>(m, "NURBS2P2D");
  add_nurbs_pyclass<2, 3>(m, "NURBS2P3D");
  add_nurbs_pyclass<2, 4>(m, "NURBS2P4D");

  add_nurbs_pyclass<3, 1>(m, "NURBS3P1D");
  add_nurbs_pyclass<3, 2>(m, "NURBS3P2D");
  add_nurbs_pyclass<3, 3>(m, "NURBS3P3D");
  add_nurbs_pyclass<3, 4>(m, "NURBS3P4D");

  add_nurbs_pyclass<4, 1>(m, "NURBS4P1D");
  add_nurbs_pyclass<4, 2>(m, "NURBS4P2D");
  add_nurbs_pyclass<4, 3>(m, "NURBS4P3D");
  add_nurbs_pyclass<4, 4>(m, "NURBS4P4D");
}
