#include <py/bezier.hpp>
#include <py/bspline.hpp>
#include <py/nurbs.hpp>

// only compile a minimal set of splines
// 3P X 3D
void init_minimal(py::module_ &m) {
  // Bezier
  add_bezier_pyclass<1, 1>(m, "Bezier1P1D");
  add_bezier_pyclass<1, 2>(m, "Bezier1P2D");
  add_bezier_pyclass<1, 3>(m, "Bezier1P3D");

  add_bezier_pyclass<2, 1>(m, "Bezier2P1D");
  add_bezier_pyclass<2, 2>(m, "Bezier2P2D");
  add_bezier_pyclass<2, 3>(m, "Bezier2P3D");

  add_bezier_pyclass<3, 1>(m, "Bezier3P1D");
  add_bezier_pyclass<3, 2>(m, "Bezier3P2D");
  add_bezier_pyclass<3, 3>(m, "Bezier3P3D");

  // BSpline
  add_bspline_pyclass<1, 1>(m, "BSpline1P1D");
  add_bspline_pyclass<1, 2>(m, "BSpline1P2D");
  add_bspline_pyclass<1, 3>(m, "BSpline1P3D");

  add_bspline_pyclass<2, 1>(m, "BSpline2P1D");
  add_bspline_pyclass<2, 2>(m, "BSpline2P2D");
  add_bspline_pyclass<2, 3>(m, "BSpline2P3D");

  add_bspline_pyclass<3, 1>(m, "BSpline3P1D");
  add_bspline_pyclass<3, 2>(m, "BSpline3P2D");
  add_bspline_pyclass<3, 3>(m, "BSpline3P3D");

  // NURBS
  add_nurbs_pyclass<1, 1>(m, "NURBS1P1D");
  add_nurbs_pyclass<1, 2>(m, "NURBS1P2D");
  add_nurbs_pyclass<1, 3>(m, "NURBS1P3D");

  add_nurbs_pyclass<2, 1>(m, "NURBS2P1D");
  add_nurbs_pyclass<2, 2>(m, "NURBS2P2D");
  add_nurbs_pyclass<2, 3>(m, "NURBS2P3D");

  add_nurbs_pyclass<3, 1>(m, "NURBS3P1D");
  add_nurbs_pyclass<3, 2>(m, "NURBS3P2D");
  add_nurbs_pyclass<3, 3>(m, "NURBS3P3D");
}
