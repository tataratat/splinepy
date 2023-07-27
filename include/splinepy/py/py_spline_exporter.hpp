#pragma once

// pybind11
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

// splinepy
#include "splinepy/py/py_spline.hpp"

/// @brief
namespace splinepy::py {

namespace py = pybind11;

/// this is a shared pointer of splinelib's SplineItem
using SplineLibIoSpline = bsplinelib::input_output::operations::SplineEntry;
// same as std::vector<SplineLibIoSpline>
using SplineLibIoSplines = bsplinelib::input_output::operations::Splines;

/// convert CoreSpline (shared_ptr of splinepybase) to splinelib io splines
/// PySpline has the same namespace
SplineLibIoSpline
PySplineToSplineLibIoSpline(std::shared_ptr<PySpline>& pyspline);

/// convert list of PySplines to vector of splinelib SplineItems
SplineLibIoSplines ListOfPySplinesToSplineLibIoSplines(py::list pysplines);

/// IGES
void ExportIges(std::string fname, py::list splines);

} // namespace splinepy::py
