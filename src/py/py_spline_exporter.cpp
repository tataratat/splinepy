#include <array>
#include <cassert>
#include <tuple>
#include <vector>

// SplineLib
#include <BSplineLib/InputOutput/iges.hpp>
#include <BSplineLib/InputOutput/operations.hpp>

// splinepy
#include "splinepy/py/py_spline_exporter.hpp"
#include "splinepy/py/py_spline_extensions.hpp"
#include "splinepy/utils/print.hpp"

/// @brief
namespace splinepy::py {

/// convert CoreSpline (shared_ptr of splinepybase) to splinelib io splines
/// PySpline has the same namespace
SplineLibIoSpline
PySplineToSplineLibIoSpline(std::shared_ptr<PySpline>& pyspline) {
  return std::dynamic_pointer_cast<typename SplineLibIoSpline::element_type>(
      splinepy::py::SameSplineWithKnotVectors(pyspline)->Core());
}

/// convert list of PySplines to vector of splinelib SplineItems
SplineLibIoSplines ListOfPySplinesToSplineLibIoSplines(py::list pysplines) {
  // prepare return obj
  SplineLibIoSplines io_splines;
  io_splines.reserve(pysplines.size());
  for (py::handle pys : pysplines) {
    io_splines.emplace_back(
        PySplineToSplineLibIoSpline(py::cast<std::shared_ptr<PySpline>&>(pys)));
  }
  return io_splines;
}

/// IGES
void ExportIges(std::string fname, py::list splines) {
  bsplinelib::input_output::iges::Write(
      ListOfPySplinesToSplineLibIoSplines(splines),
      fname);
}

/// @brief Adds Python spline exporter
/// @param m Python module
void init_spline_exporter(py::module& m) {

  m.def("export_iges",
        &splinepy::py::ExportIges,
        py::arg("fname"),
        py::arg("splines"));
}

} // namespace splinepy::py
