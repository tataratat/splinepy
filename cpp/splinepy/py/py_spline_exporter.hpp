#pragma once

#include <array>
#include <cassert>
#include <tuple>
#include <vector>

// pybind11
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

// SplineLib
#include <Sources/InputOutput/iges.hpp>
#include <Sources/InputOutput/irit.hpp>
#include <Sources/InputOutput/operations.hpp>
#include <Sources/InputOutput/vtk.hpp>
#include <Sources/InputOutput/xml.hpp>

// splinepy
#include <splinepy/py/py_spline.hpp>
#include <splinepy/py/py_spline_extensions.hpp>
#include <splinepy/utils/print.hpp>

/// @brief
namespace splinepy::py {

namespace py = pybind11;

/// this is a shared pointer of splinelib's SplineItem
using SplineLibIoSpline =
    splinelib::sources::input_output::operations::SplineEntry;
// same as std::vector<SplineLibIoSpline>
using SplineLibIoSplines =
    splinelib::sources::input_output::operations::Splines;

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
  splinelib::sources::input_output::iges::Write(
      ListOfPySplinesToSplineLibIoSplines(splines),
      fname);
}

/// IRIT
void ExportIrit(std::string fname, py::list splines) {
  splinelib::sources::input_output::irit::Write(
      ListOfPySplinesToSplineLibIoSplines(splines),
      fname);
}

/// XML
void ExportXml(std::string fname, py::list splines) {
  splinelib::sources::input_output::xml::Write(
      ListOfPySplinesToSplineLibIoSplines(splines),
      fname);
}

/// VTK - sampled spline
void ExportVtk(std::string fname,
               py::list splines,
               py::list resolutions_per_spline) {
  using ResolutionsPerSpline =
      splinelib::sources::input_output::vtk::NumbersOfParametricCoordinates;
  using ResolutionsType = typename ResolutionsPerSpline::value_type;
  using ResolutionValueType =
      typename ResolutionsPerSpline::value_type::value_type;

  // quick check
  if (splines.size() != resolutions_per_spline.size()) {
    splinepy::utils::PrintAndThrowError("Number of splines (",
                                        splines.size(),
                                        ") and number of resolutions (",
                                        resolutions_per_spline.size(),
                                        ") does not match.");
  }

  // get vector of splinelib io splines
  auto sl_io_splines = ListOfPySplinesToSplineLibIoSplines(splines);

  // prepare vector of resolutions
  ResolutionsPerSpline sl_rps;
  sl_rps.reserve(resolutions_per_spline.size());

  // loop resolutions
  int i{};
  for (py::handle res : resolutions_per_spline) {
    py::array_t<int> r = py::cast<py::array_t<int>>(res);
    int* r_ptr = static_cast<int*>(r.request().ptr);
    const int r_len = r.size();

    // check resolution len
    if (r_len != sl_io_splines[i]->parametric_dimensionality_) {
      splinepy::utils::PrintAndThrowError(
          "Invalid resolutions length for spline index (",
          i,
          ").",
          "Expected: ",
          sl_io_splines[i]->parametric_dimensionality_,
          "Given: ",
          r_len);
    }

    // prepare resolutions
    ResolutionsType sl_res;
    sl_res.reserve(r_len);
    for (int j{}; j < r_len; ++j) {
      // check if values are at least 2
      const int& res_value = r_ptr[j];
      if (res_value < 2) {
        splinepy::utils::PrintAndThrowError(
            "Invalid resolution value at index (",
            j,
            ")",
            "for spline index (",
            i,
            ").",
            "Expected to be greater than 2. Given:",
            res_value);
      }
      sl_res.emplace_back(ResolutionValueType{res_value});
    }
    // append resolutions
    sl_rps.push_back(std::move(sl_res));

    // prepare next loop
    ++i;
  }

  // good,
  splinelib::sources::input_output::vtk::Sample(sl_io_splines, fname, sl_rps);
}

/// @brief Adds Python spline exporter
/// @param m Python module
inline void add_spline_exporter(py::module& m) {

  m.def("export_iges",
        &splinepy::py::ExportIges,
        py::arg("fname"),
        py::arg("splines"));
  m.def("export_irit",
        &splinepy::py::ExportIrit,
        py::arg("fname"),
        py::arg("splines"));
  m.def("export_xml",
        &splinepy::py::ExportXml,
        py::arg("fname"),
        py::arg("splines"));
  m.def("export_vtk",
        &splinepy::py::ExportVtk,
        py::arg("fname"),
        py::arg("splines"),
        py::arg("resolutions_per_spline"));
}

} // namespace splinepy::py
