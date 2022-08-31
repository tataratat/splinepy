#include <pybind11/pybind11.h>

namespace py = pybind11;

// Bezier
void init_bezier1(py::module_&);
void init_bezier2(py::module_&);
void init_bezier3(py::module_&);
void init_bezier4(py::module_&);
void init_bezier5(py::module_&);
void init_bezier6(py::module_&);
void init_bezier7(py::module_&);
void init_bezier8(py::module_&);
void init_bezier9(py::module_&);
void init_bezier10(py::module_&);

// Rational Bezier
void init_rational_bezier1(py::module_&);
void init_rational_bezier2(py::module_&);
void init_rational_bezier3(py::module_&);
void init_rational_bezier4(py::module_&);
void init_rational_bezier5(py::module_&);
void init_rational_bezier6(py::module_&);
void init_rational_bezier7(py::module_&);
void init_rational_bezier8(py::module_&);
void init_rational_bezier9(py::module_&);
void init_rational_bezier10(py::module_&);

// BSpline
void init_bspline1(py::module_&);
void init_bspline2(py::module_&);
void init_bspline3(py::module_&);
void init_bspline4(py::module_&);
void init_bspline5(py::module_&);
void init_bspline6(py::module_&);
void init_bspline7(py::module_&);
void init_bspline8(py::module_&);
void init_bspline9(py::module_&);
void init_bspline10(py::module_&);

// NURBS
void init_nurbs1(py::module_&);
void init_nurbs2(py::module_&);
void init_nurbs3(py::module_&);
void init_nurbs4(py::module_&);
void init_nurbs5(py::module_&);
void init_nurbs6(py::module_&);
void init_nurbs7(py::module_&);
void init_nurbs8(py::module_&);
void init_nurbs9(py::module_&);
void init_nurbs10(py::module_&);

// minimal
void init_minimal(py::module_&);

// Reader
void init_reader(py::module_&);

// Exporter
void init_exporter(py::module_&);

PYBIND11_MODULE(_splinepy, m) {
#ifdef _MINIMAL_
  init_minimal(m);
#else
  init_bezier1(m);
  init_bezier2(m);
  init_bezier3(m);
  init_bezier4(m);
  init_bezier5(m);
  init_bezier6(m);
  init_bezier7(m);
  init_bezier8(m);
  init_bezier9(m);
  init_bezier10(m);

  init_rational_bezier1(m);
  init_rational_bezier2(m);
  init_rational_bezier3(m);
  init_rational_bezier4(m);
  init_rational_bezier5(m);
  init_rational_bezier6(m);
  init_rational_bezier7(m);
  init_rational_bezier8(m);
  init_rational_bezier9(m);
  init_rational_bezier10(m);

  init_bspline1(m);
  init_bspline2(m);
  init_bspline3(m);
  init_bspline4(m);
  init_bspline5(m);
  init_bspline6(m);
  init_bspline7(m);
  init_bspline8(m);
  init_bspline9(m);
  init_bspline10(m);

  init_nurbs1(m);
  init_nurbs2(m);
  init_nurbs3(m);
  init_nurbs4(m);
  init_nurbs5(m);
  init_nurbs6(m);
  init_nurbs7(m);
  init_nurbs8(m);
  init_nurbs9(m);
  init_nurbs10(m);
#endif

  init_reader(m);
  init_exporter(m);
}
