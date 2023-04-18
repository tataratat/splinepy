#pragma once

#include <array>
#include <cassert>
#include <tuple>
#include <vector>

// pybind11
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

// Bezman
#include <bezman/src/utils/algorithms/point_uniquifier.hpp>

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

template<std::size_t parametric_dimension, std::size_t physical_dimension>
py::array_t<int>
InterfacesFromBoundaryCenters_(const py::array_t<double>& py_center_vertices,
                               const double& tolerance) {
  // Auxiliary Function to reduce total number of declarations
  using PhysicalPointType = bezman::Point<physical_dimension, double>;

  // Determine data
  double* centers_ptr = static_cast<double*>(py_center_vertices.request().ptr);
  const std::size_t number_of_center_points =
      py_center_vertices.request().shape[0];
  const std::size_t physical_dimension_ = py_center_vertices.request().shape[1];
  constexpr std::size_t n_faces_per_patch = parametric_dimension * 2;

  // Assertions
  assert(number_of_center_points > 0);
  assert(physical_dimension_ > 0);
  assert(number_of_center_points % n_faces_per_patch == 0);

  // Convert points into bezman points
  std::vector<PhysicalPointType> center_points;

  PhysicalPointType minimumVertex{}, maximumVertex{};
  // Assign first vertex to both min and max
  for (std::size_t i_dim{}; i_dim < physical_dimension; i_dim++) {
    minimumVertex[i_dim] = centers_ptr[i_dim];
    maximumVertex[i_dim] = centers_ptr[i_dim];
  }

  center_points.reserve(number_of_center_points);
  for (std::size_t i_point{}; i_point < number_of_center_points; i_point++) {
    PhysicalPointType point{};
    for (std::size_t i_dim{}; i_dim < physical_dimension; i_dim++) {
      point[i_dim] = centers_ptr[i_point * physical_dimension_ + i_dim];
      minimumVertex[i_dim] =
          std::min(minimumVertex[i_dim],
                   centers_ptr[i_point * physical_dimension_ + i_dim]);
      maximumVertex[i_dim] =
          std::max(maximumVertex[i_dim],
                   centers_ptr[i_point * physical_dimension_ + i_dim]);
    }
    center_points.push_back(point);
  }

  // Hand to bezman for connectivity
  const auto connectivity =
      bezman::utils::algorithms::FindConnectivityFromCenters<
          parametric_dimension,
          false>(center_points, maximumVertex - minimumVertex, tolerance);

  // Transform points into an array
  const int number_of_patches = connectivity.size();
  py::array_t<int> py_connectivity =
      py::array_t<int>(number_of_patches * n_faces_per_patch);
  py_connectivity.resize({(int) number_of_patches, (int) n_faces_per_patch});
  int* py_connectivity_ptr = static_cast<int*>(py_connectivity.request().ptr);
  for (std::size_t i_patch{}; i_patch < connectivity.size(); i_patch++) {
    for (std::size_t i_face{}; i_face < n_faces_per_patch; i_face++) {
      py_connectivity_ptr[i_patch * n_faces_per_patch + i_face] =
          static_cast<int>(connectivity[i_patch][i_face]);
    }
  }

  return py_connectivity;
}

/**
 * @brief  Determines the Connectivity of spline patches
 *
 * @param py_center_vertices  Vertices in the center of the boundaries
 * @param tolerance tolerance between two neighboring face centers for them
 * to be fused
 * @param parametric_dimension Parametric dimension of the spline grid
 * @return py::array_t<int> connectivity
 */
py::array_t<int>
InterfacesFromBoundaryCenters(const py::array_t<double>& py_center_vertices,
                              const double& tolerance,
                              const int& parametric_dimension) {
  // Transform points from pyarray into bezman point vector
  double* centers_ptr = static_cast<double*>(py_center_vertices.request().ptr);
  const std::size_t physical_dimension_ = py_center_vertices.request().shape[1];
  const std::size_t number_of_center_points =
      py_center_vertices.request().shape[0];

  // Check input data
  assert(0 == (number_of_center_points % (2 * parametric_dimension)));

  // Convert points into bezman type points
  switch (physical_dimension_) {
  case 1:
    switch (parametric_dimension) {
    case 1:
      return InterfacesFromBoundaryCenters_<1uL, 1uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 2:
      return InterfacesFromBoundaryCenters_<2uL, 1uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 3:
      return InterfacesFromBoundaryCenters_<3uL, 1uL>(py_center_vertices,
                                                      tolerance);
      break;
#ifdef SPLINEPY_MORE
    case 4:
      return InterfacesFromBoundaryCenters_<4uL, 1uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 5:
      return InterfacesFromBoundaryCenters_<5uL, 1uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 6:
      return InterfacesFromBoundaryCenters_<6uL, 1uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 7:
      return InterfacesFromBoundaryCenters_<7uL, 1uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 8:
      return InterfacesFromBoundaryCenters_<8uL, 1uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 9:
      return InterfacesFromBoundaryCenters_<9uL, 1uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 10:
      return InterfacesFromBoundaryCenters_<10uL, 1uL>(py_center_vertices,
                                                       tolerance);
      break;
#endif
    default:
      break;
    }
    break;
  case 2:
    switch (parametric_dimension) {
    case 1:
      return InterfacesFromBoundaryCenters_<1uL, 2uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 2:
      return InterfacesFromBoundaryCenters_<2uL, 2uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 3:
      return InterfacesFromBoundaryCenters_<3uL, 2uL>(py_center_vertices,
                                                      tolerance);
      break;
#ifdef SPLINEPY_MORE
    case 4:
      return InterfacesFromBoundaryCenters_<4uL, 2uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 5:
      return InterfacesFromBoundaryCenters_<5uL, 2uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 6:
      return InterfacesFromBoundaryCenters_<6uL, 2uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 7:
      return InterfacesFromBoundaryCenters_<7uL, 2uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 8:
      return InterfacesFromBoundaryCenters_<8uL, 2uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 9:
      return InterfacesFromBoundaryCenters_<9uL, 2uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 10:
      return InterfacesFromBoundaryCenters_<10uL, 2uL>(py_center_vertices,
                                                       tolerance);
      break;
#endif

    default:
      break;
    }
    break;
  case 3:
    switch (parametric_dimension) {
    case 1:
      return InterfacesFromBoundaryCenters_<1uL, 3uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 2:
      return InterfacesFromBoundaryCenters_<2uL, 3uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 3:
      return InterfacesFromBoundaryCenters_<3uL, 3uL>(py_center_vertices,
                                                      tolerance);
      break;
#ifdef SPLINEPY_MORE
    case 4:
      return InterfacesFromBoundaryCenters_<4uL, 3uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 5:
      return InterfacesFromBoundaryCenters_<5uL, 3uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 6:
      return InterfacesFromBoundaryCenters_<6uL, 3uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 7:
      return InterfacesFromBoundaryCenters_<7uL, 3uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 8:
      return InterfacesFromBoundaryCenters_<8uL, 3uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 9:
      return InterfacesFromBoundaryCenters_<9uL, 3uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 10:
      return InterfacesFromBoundaryCenters_<10uL, 3uL>(py_center_vertices,
                                                       tolerance);
      break;
#endif

    default:
      break;
    }
    break;
#ifdef SPLINEPY_MORE
  case 4:
    switch (parametric_dimension) {
    case 1:
      return InterfacesFromBoundaryCenters_<1uL, 4uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 2:
      return InterfacesFromBoundaryCenters_<2uL, 4uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 3:
      return InterfacesFromBoundaryCenters_<3uL, 4uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 4:
      return InterfacesFromBoundaryCenters_<4uL, 4uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 5:
      return InterfacesFromBoundaryCenters_<5uL, 4uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 6:
      return InterfacesFromBoundaryCenters_<6uL, 4uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 7:
      return InterfacesFromBoundaryCenters_<7uL, 4uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 8:
      return InterfacesFromBoundaryCenters_<8uL, 4uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 9:
      return InterfacesFromBoundaryCenters_<9uL, 4uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 10:
      return InterfacesFromBoundaryCenters_<10uL, 4uL>(py_center_vertices,
                                                       tolerance);
      break;
    default:
      break;
    }
  case 5:
    switch (parametric_dimension) {
    case 1:
      return InterfacesFromBoundaryCenters_<1uL, 5uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 2:
      return InterfacesFromBoundaryCenters_<2uL, 5uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 3:
      return InterfacesFromBoundaryCenters_<3uL, 5uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 4:
      return InterfacesFromBoundaryCenters_<4uL, 5uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 5:
      return InterfacesFromBoundaryCenters_<5uL, 5uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 6:
      return InterfacesFromBoundaryCenters_<6uL, 5uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 7:
      return InterfacesFromBoundaryCenters_<7uL, 5uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 8:
      return InterfacesFromBoundaryCenters_<8uL, 5uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 9:
      return InterfacesFromBoundaryCenters_<9uL, 5uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 10:
      return InterfacesFromBoundaryCenters_<10uL, 5uL>(py_center_vertices,
                                                       tolerance);
      break;
    default:
      break;
    }
  case 6:
    switch (parametric_dimension) {
    case 1:
      return InterfacesFromBoundaryCenters_<1uL, 6uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 2:
      return InterfacesFromBoundaryCenters_<2uL, 6uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 3:
      return InterfacesFromBoundaryCenters_<3uL, 6uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 4:
      return InterfacesFromBoundaryCenters_<4uL, 6uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 5:
      return InterfacesFromBoundaryCenters_<5uL, 6uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 6:
      return InterfacesFromBoundaryCenters_<6uL, 6uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 7:
      return InterfacesFromBoundaryCenters_<7uL, 6uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 8:
      return InterfacesFromBoundaryCenters_<8uL, 6uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 9:
      return InterfacesFromBoundaryCenters_<9uL, 6uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 10:
      return InterfacesFromBoundaryCenters_<10uL, 6uL>(py_center_vertices,
                                                       tolerance);
      break;
    default:
      break;
    }
  case 7:
    switch (parametric_dimension) {
    case 1:
      return InterfacesFromBoundaryCenters_<1uL, 7uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 2:
      return InterfacesFromBoundaryCenters_<2uL, 7uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 3:
      return InterfacesFromBoundaryCenters_<3uL, 7uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 4:
      return InterfacesFromBoundaryCenters_<4uL, 7uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 5:
      return InterfacesFromBoundaryCenters_<5uL, 7uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 6:
      return InterfacesFromBoundaryCenters_<6uL, 7uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 7:
      return InterfacesFromBoundaryCenters_<7uL, 7uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 8:
      return InterfacesFromBoundaryCenters_<8uL, 7uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 9:
      return InterfacesFromBoundaryCenters_<9uL, 7uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 10:
      return InterfacesFromBoundaryCenters_<10uL, 7uL>(py_center_vertices,
                                                       tolerance);
      break;
    default:
      break;
    }
  case 8:
    switch (parametric_dimension) {
    case 1:
      return InterfacesFromBoundaryCenters_<1uL, 8uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 2:
      return InterfacesFromBoundaryCenters_<2uL, 8uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 3:
      return InterfacesFromBoundaryCenters_<3uL, 8uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 4:
      return InterfacesFromBoundaryCenters_<4uL, 8uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 5:
      return InterfacesFromBoundaryCenters_<5uL, 8uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 6:
      return InterfacesFromBoundaryCenters_<6uL, 8uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 7:
      return InterfacesFromBoundaryCenters_<7uL, 8uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 8:
      return InterfacesFromBoundaryCenters_<8uL, 8uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 9:
      return InterfacesFromBoundaryCenters_<9uL, 8uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 10:
      return InterfacesFromBoundaryCenters_<10uL, 8uL>(py_center_vertices,
                                                       tolerance);
      break;
    default:
      break;
    }
  case 9:
    switch (parametric_dimension) {
    case 1:
      return InterfacesFromBoundaryCenters_<1uL, 9uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 2:
      return InterfacesFromBoundaryCenters_<2uL, 9uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 3:
      return InterfacesFromBoundaryCenters_<3uL, 9uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 4:
      return InterfacesFromBoundaryCenters_<4uL, 9uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 5:
      return InterfacesFromBoundaryCenters_<5uL, 9uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 6:
      return InterfacesFromBoundaryCenters_<6uL, 9uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 7:
      return InterfacesFromBoundaryCenters_<7uL, 9uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 8:
      return InterfacesFromBoundaryCenters_<8uL, 9uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 9:
      return InterfacesFromBoundaryCenters_<9uL, 9uL>(py_center_vertices,
                                                      tolerance);
      break;
    case 10:
      return InterfacesFromBoundaryCenters_<10uL, 9uL>(py_center_vertices,
                                                       tolerance);
      break;
    default:
      break;
    }
  case 10:
    switch (parametric_dimension) {
    case 1:
      return InterfacesFromBoundaryCenters_<1uL, 10uL>(py_center_vertices,
                                                       tolerance);
      break;
    case 2:
      return InterfacesFromBoundaryCenters_<2uL, 10uL>(py_center_vertices,
                                                       tolerance);
      break;
    case 3:
      return InterfacesFromBoundaryCenters_<3uL, 10uL>(py_center_vertices,
                                                       tolerance);
      break;
    case 4:
      return InterfacesFromBoundaryCenters_<4uL, 10uL>(py_center_vertices,
                                                       tolerance);
      break;
    case 5:
      return InterfacesFromBoundaryCenters_<5uL, 10uL>(py_center_vertices,
                                                       tolerance);
      break;
    case 6:
      return InterfacesFromBoundaryCenters_<6uL, 10uL>(py_center_vertices,
                                                       tolerance);
      break;
    case 7:
      return InterfacesFromBoundaryCenters_<7uL, 10uL>(py_center_vertices,
                                                       tolerance);
      break;
    case 8:
      return InterfacesFromBoundaryCenters_<8uL, 10uL>(py_center_vertices,
                                                       tolerance);
      break;
    case 9:
      return InterfacesFromBoundaryCenters_<9uL, 10uL>(py_center_vertices,
                                                       tolerance);
      break;
    case 10:
      return InterfacesFromBoundaryCenters_<10uL, 10uL>(py_center_vertices,
                                                        tolerance);
      break;
    default:
      break;
    }
#endif
  default:
    break;
  }

#ifdef SPLINEPY_MORE
  splinepy::utils::PrintAndThrowError(
      "Only implemented for <2-10> : <2-10> dimensions");
#else
  splinepy::utils::PrintAndThrowError(
      "Only implemented for <1-3> : <1-3> dimensions");
#endif
  // dummy statement for compiler
  return py::array_t<int>();
}

/**
 * @brief Orientation between two adjacent splines
 *
 * If two splines share the same boundary this function retrieves their
 * orientation, by mapping the mappings of the parametric axis onto each other.
 * This is (among others) required for Gismo and Nutils export
 *
 * @param pyspline_start Spline object from start
 * @param boundary_start Boundary ID from start spline
 * @param pyspline_end Spline object from end *to which is mapped
 * @param boundary_end Boundary ID of adjacent spline
 * @param int_mappings_ptr (output) integer mappings
 * @param bool_orientations_ptr (output) axis alignement
 * @return void
 */
void GetBoundaryOrientation(
    const std::shared_ptr<splinepy::splines::SplinepyBase>& pyspline_start,
    const int& boundary_start,
    const std::shared_ptr<splinepy::splines::SplinepyBase>& pyspline_end,
    const int& boundary_end,
    const double& tolerance,
    int* int_mappings_ptr,
    bool* bool_orientations_ptr) {
  // Init return values and get auxiliary data
  const int& para_dim_ = pyspline_start->SplinepyParaDim();
  const int& dim_ = pyspline_start->SplinepyDim();

  // Checks
  if ((para_dim_ != pyspline_end->SplinepyParaDim())
      || (dim_ != pyspline_end->SplinepyDim())) {
    splinepy::utils::PrintAndThrowError(
        "Spline Orientation can not be checked, as they have mismatching"
        "dimensionality start spline has dimensions ",
        para_dim_,
        "D -> ",
        dim_,
        "D, the adjacent one has dimensions ",
        pyspline_end->SplinepyParaDim(),
        "D -> ",
        pyspline_end->SplinepyDim(),
        "D.");
  }

  // First Check the orientation of the first entry by comparing their ids
  const int boundary_start_p_dim = static_cast<int>(boundary_start / 2);
  const bool boundary_start_orientation = (boundary_start % 2) == 0;
  const int boundary_end_p_dim = static_cast<int>(boundary_end / 2);
  const bool boundary_end_orientation = (boundary_end % 2) == 0;
  int_mappings_ptr[boundary_start_p_dim] = boundary_end_p_dim;
  // Note: Here might be a discrepency with gismo's orientation, and it needs to
  // be checked in the future. I am awaiting a response from gismo developers,
  // it is poosible the orientation of the interface edge might be flipped
  // (bugfix: negate the following expression)
  bool_orientations_ptr[boundary_start_p_dim] =
      (boundary_start_orientation ^ boundary_end_orientation);

  /// Compare jacobians for remaining entries
  // Calculate Parametric bounds
  std::vector<double> bounds_start(para_dim_ * 2);
  pyspline_start->SplinepyParametricBounds(bounds_start.data());
  std::vector<double> bounds_end(para_dim_ * 2);
  pyspline_end->SplinepyParametricBounds(bounds_end.data());
  // Determine face center position in parametric space
  std::vector<double> boundary_center_start(para_dim_),
      boundary_center_end(para_dim_);

  for (int i{}; i < para_dim_; i++) {
    if (i == boundary_start_p_dim) {
      boundary_center_start[i] = boundary_start_orientation
                                     ? bounds_start[i]
                                     : bounds_start[i + para_dim_];
    } else {
      boundary_center_start[i] =
          .5 * (bounds_start[i] + bounds_start[i + para_dim_]);
    }
    if (i == boundary_end_p_dim) {
      boundary_center_end[i] =
          boundary_end_orientation ? bounds_end[i] : bounds_end[i + para_dim_];
    } else {
      boundary_center_end[i] = .5 * (bounds_end[i] + bounds_end[i + para_dim_]);
    }
  }

  // Calculate Jacobians
  std::vector<double> jacobian_start(para_dim_ * dim_),
      jacobian_end(para_dim_ * dim_);
  pyspline_start->SplinepyJacobian(boundary_center_start.data(),
                                   jacobian_start.data());
  pyspline_end->SplinepyJacobian(boundary_center_end.data(),
                                 jacobian_end.data());

  // Check the angle between the jacobian entries
  for (int i_pd{}; i_pd < para_dim_; i_pd++) {
    if (i_pd == boundary_start_p_dim) {
      continue;
    }
    double norm_s{};
    for (int k{}; k < dim_; k++) {
      // [i_query * pdim * dim + i_paradim * dim + i_dim]
      norm_s += jacobian_start[i_pd + k * para_dim_]
                * jacobian_start[i_pd + k * para_dim_];
    }
    for (int j{}; j < para_dim_; j++) {
      double norm_e{}, dot_p{};
      for (int k{}; k < dim_; k++) {
        dot_p += jacobian_start[i_pd + k * para_dim_]
                 * jacobian_end[j + k * para_dim_];
        norm_e +=
            jacobian_end[j + k * para_dim_] * jacobian_end[j + k * para_dim_];
      }

      // Check angle
      const double cos_angle = abs(dot_p / std::sqrt(norm_s * norm_e));
      if (cos_angle > (1. - tolerance)) {
        int_mappings_ptr[i_pd] = j;
        bool_orientations_ptr[i_pd] = (dot_p > 0);
        break;
      }
    }
  }
}

/**
 * @brief Get the Boundary Orientations object
 *
 * @param spline_list
 * @param base_id
 * @param base_face_id
 * @param base_id
 * @param base_face_id
 * @param tolerance
 * @param n_threads
 * @return py::tuple
 */
py::tuple GetBoundaryOrientations(const py::list& spline_list,
                                  const py::array_t<int>& base_id,
                                  const py::array_t<int>& base_face_id,
                                  const py::array_t<int>& neighbor_id,
                                  const py::array_t<int>& neighbor_face_id,
                                  const double tolerance,
                                  const int n_threads) {
  // Basic Checks
  // Check if all have same size
  if (!((base_id.size() == base_face_id.size())
        && (neighbor_id.size() == neighbor_face_id.size())
        && (base_id.size() == neighbor_face_id.size()))) {
    splinepy::utils::PrintAndThrowError(
        "The ID arrays need to be of same size, please check for "
        "consistencies.");
  }

  // Auxiliary data
  const int* base_id_ptr = static_cast<int*>(base_id.request().ptr);
  const int* base_face_id_ptr = static_cast<int*>(base_face_id.request().ptr);
  const int* neighbor_id_ptr = static_cast<int*>(neighbor_id.request().ptr);
  const int* neighbor_face_id_ptr =
      static_cast<int*>(neighbor_face_id.request().ptr);
  const auto cpp_spline_list =
      ListOfPySplinesToVectorOfCoreSplines(spline_list);
  const int n_connections = base_id.size();

  const int para_dim_ = cpp_spline_list[0]->SplinepyParaDim();

  py::array_t<int> int_mapping(n_connections * para_dim_);
  int* int_mapping_ptr = static_cast<int*>(int_mapping.request().ptr);
  py::array_t<bool> bool_orientations(n_connections * para_dim_);
  bool* bool_orientations_ptr =
      static_cast<bool*>(bool_orientations.request().ptr);

  // Provide lambda for multithread execution
  auto get_orientation = [&](int start, int end) {
    for (int i{start}; i < end; ++i) {
      GetBoundaryOrientation(cpp_spline_list[base_id_ptr[i]],
                             base_face_id_ptr[i],
                             cpp_spline_list[neighbor_id_ptr[i]],
                             neighbor_face_id_ptr[i],
                             tolerance,
                             &int_mapping_ptr[i * para_dim_],
                             &bool_orientations_ptr[i * para_dim_]);
    }
  };

  // Execute in parallel
  splinepy::utils::NThreadExecution(get_orientation, n_connections, n_threads);

  // Resize and return
  int_mapping.resize({n_connections, para_dim_});
  bool_orientations.resize({n_connections, para_dim_});

  return py::make_tuple(int_mapping, bool_orientations);
}

/**
 * @brief Retrieve information related to mfem export
 *
 * @param py_corner_vertices vertices at the spline-corners
 * @param tolerance tolerance to delete duplicates
 * @return py::tuple with
 *    py::array_t<int> : connectivity
 *    py::array_t<int> : vertex_ids
 *    py::array_t<int> : edges
 *    py::array_t<int> : boundaries
 *    bool             : is structured mesh
 */
py::tuple RetrieveMfemInformation(const py::array_t<double>& py_corner_vertices,
                                  const double& tolerance) {
  // Unfortunatly bezman requires point-types to perform routines and does not
  // work on py arrays All of the arguments serve as outputs except for
  // corner_vertices
  py::buffer_info corner_vertex_buffer = py_corner_vertices.request();
  double* corner_ptr = static_cast<double*>(corner_vertex_buffer.ptr);
  const std::size_t physical_dimension_ = corner_vertex_buffer.shape[1];
  const std::size_t number_of_corner_points = corner_vertex_buffer.shape[0];
  // Check if mesh can be used for mfem mesh
  bool is_structured{true};
  if (physical_dimension_ == 2) {
    // Check if vertex size checks out
    assert(number_of_corner_points % 4 == 0);
    const std::size_t number_of_patches = number_of_corner_points / 4;

    // Transform corner_vertex into arrays
    using PointType = bezman::Point<2ul, double>;
    std::vector<PointType> corner_vertices;
    corner_vertices.reserve(number_of_corner_points);
    // Init Max and min vertices to determine metric
    PointType maxvertex{corner_ptr[0], corner_ptr[1]},
        minvertex{corner_ptr[0], corner_ptr[1]};
    for (std::size_t i_c{}; i_c < number_of_corner_points; i_c++) {
      PointType vertex{};
      for (std::size_t i_dim{}; i_dim < physical_dimension_; i_dim++) {
        vertex[i_dim] = corner_ptr[i_c * physical_dimension_ + i_dim];
        maxvertex[i_dim] = std::max(vertex[i_dim], maxvertex[i_dim]);
        minvertex[i_dim] = std::min(vertex[i_dim], minvertex[i_dim]);
      }
      corner_vertices.push_back(vertex);
    }

    // Retrieve MFEM information using bezman
    // Connectivity :     std::vector<std::array<std::size_t, 4>>
    // Vertex_ids :       std::vector<std::size_t>
    // edge_information : std::vector<std::array<std::size_t,3>>
    // boundaries :       std::vector<std::array<std::size_t,2>>
    const auto [connectivity, vertex_ids, edge_information, boundaries] =
        [&]() {
          try {
            return bezman::utils::algorithms::ExtractMFEMInformation(
                corner_vertices,
                maxvertex - minvertex,
                tolerance);
          } catch (...) {
            is_structured = false;
            return std::make_tuple(
                // Connectivity
                bezman::utils::algorithms::FindConnectivityFromCorners<2>(
                    corner_vertices,
                    maxvertex - minvertex,
                    tolerance,
                    false),
                // All others initialized empty
                std::vector<std::size_t>{},
                std::vector<std::array<std::size_t, 3>>{},
                std::vector<std::array<std::size_t, 2>>{});
          }
        }();

    // -- Transform data to python format --
    // Connectivity
    assert(connectivity.size() == number_of_patches);
    py::array_t<int> py_connectivity =
        py::array_t<int>(connectivity.size() * 4);
    py_connectivity.resize({(int) number_of_patches, (int) 4});
    int* py_connectivity_ptr = static_cast<int*>(py_connectivity.request().ptr);
    for (std::size_t i_patch{}; i_patch < connectivity.size(); i_patch++) {
      for (std::size_t i_face{}; i_face < 4ul; i_face++) {
        py_connectivity_ptr[i_patch * 4ul + i_face] =
            static_cast<int>(connectivity[i_patch][i_face]);
      }
    }

    // Return only connectivity if the mesh is unstructured and can not be used
    // for mfem export
    if (!is_structured) {
      return py::make_tuple(py_connectivity,
                            py::array_t<int>{},
                            py::array_t<int>{},
                            py::array_t<int>{},
                            is_structured);
    }

    // Vertex IDS
    assert(vertex_ids.size() == corner_vertices.size());
    py::array_t<int> py_vertex_ids = py::array_t<int>(number_of_corner_points);
    py_vertex_ids.resize({(int) number_of_corner_points});
    int* py_vertex_ids_ptr = static_cast<int*>(py_vertex_ids.request().ptr);
    for (std::size_t i_ctps{}; i_ctps < number_of_corner_points; i_ctps++) {
      py_vertex_ids_ptr[i_ctps] = static_cast<int>(vertex_ids[i_ctps]);
    }

    // Edges
    assert(edge_information.size() > 0);
    py::array_t<int> py_edges = py::array_t<int>(edge_information.size() * 3);
    py_edges.resize({(int) edge_information.size(), (int) 3});
    int* py_edges_ptr = static_cast<int*>(py_edges.request().ptr);
    for (std::size_t i_edge{}; i_edge < edge_information.size(); i_edge++) {
      for (std::size_t i_face{}; i_face < 3ul; i_face++) {
        py_edges_ptr[i_edge * 3ul + i_face] =
            static_cast<int>(edge_information[i_edge][i_face]);
      }
    }

    // Boundaries
    assert(boundaries.size() > 0);
    py::array_t<int> py_boundaries = py::array_t<int>(boundaries.size() * 2);
    py_boundaries.resize({(int) boundaries.size(), (int) 2});
    int* py_boundaries_ptr = static_cast<int*>(py_boundaries.request().ptr);
    for (std::size_t i_boundary{}; i_boundary < boundaries.size();
         i_boundary++) {
      for (std::size_t i_id{}; i_id < 2ul; i_id++) {
        py_boundaries_ptr[i_boundary * 2ul + i_id] =
            static_cast<int>(boundaries[i_boundary][i_id]);
      }
    }

    return py::make_tuple(py_connectivity,
                          py_vertex_ids,
                          py_edges,
                          py_boundaries,
                          is_structured);

  } else if (physical_dimension_ == 3) {
    // Check if vertex size checks out
    assert(number_of_corner_points % 8 == 0);
    const std::size_t number_of_patches = number_of_corner_points / 8;

    // Transform corner_vertex into arrays
    using PointType = bezman::Point<3ul, double>;
    std::vector<PointType> corner_vertices;
    corner_vertices.reserve(number_of_corner_points);
    // Init Max and min vertices to determine metric
    PointType maxvertex{corner_ptr[0], corner_ptr[1], corner_ptr[2]},
        minvertex{corner_ptr[0], corner_ptr[1], corner_ptr[2]};
    for (std::size_t i_c{}; i_c < number_of_corner_points; i_c++) {
      PointType vertex{};
      for (std::size_t i_dim{}; i_dim < physical_dimension_; i_dim++) {
        vertex[i_dim] = corner_ptr[i_c * physical_dimension_ + i_dim];
        maxvertex[i_dim] = std::max(vertex[i_dim], maxvertex[i_dim]);
        minvertex[i_dim] = std::min(vertex[i_dim], minvertex[i_dim]);
      }
      corner_vertices.push_back(vertex);
    }
    // Retrieve MFEM information using bezman
    // Connectivity : std::vector<std::array<std::size_t, 6>>
    // Vertex_ids : std::vector<std::std::size_t>
    // edge_information : std::vector<std::array<std::size_t,3>
    // boundaries : std::vector<std::array<std::size_t,4>
    const auto [connectivity, vertex_ids, edge_information, boundaries] =
        [&]() {
          try {
            return bezman::utils::algorithms::ExtractMFEMInformation(
                corner_vertices,
                maxvertex - minvertex,
                tolerance);
          } catch (...) {
            is_structured = false;
            return std::make_tuple(
                // Connectivity
                bezman::utils::algorithms::FindConnectivityFromCorners<3>(
                    corner_vertices,
                    maxvertex - minvertex,
                    tolerance,
                    false),
                // All others initialized empty
                std::vector<std::size_t>{},
                std::vector<std::array<std::size_t, 3>>{},
                std::vector<std::array<std::size_t, 4>>{});
          }
        }();

    // -- Transform data to python format --
    // Connectivity
    assert(connectivity.size() == number_of_patches);
    py::array_t<int> py_connectivity =
        py::array_t<int>(connectivity.size() * 6);
    py_connectivity.resize({(int) number_of_patches, (int) 6});
    int* py_connectivity_ptr = static_cast<int*>(py_connectivity.request().ptr);
    for (std::size_t i_patch{}; i_patch < connectivity.size(); i_patch++) {
      for (std::size_t i_face{}; i_face < 6ul; i_face++) {
        py_connectivity_ptr[i_patch * 6ul + i_face] =
            static_cast<int>(connectivity[i_patch][i_face]);
      }
    }

    // Return only connectivity if the mesh is unstructured and can not be used
    // for mfem export
    if (!is_structured) {
      return py::make_tuple(py_connectivity,
                            py::array_t<int>{},
                            py::array_t<int>{},
                            py::array_t<int>{},
                            is_structured);
    }

    // Vertex IDS
    assert(vertex_ids.size() == corner_vertices.size());
    py::array_t<int> py_vertex_ids = py::array_t<int>(number_of_corner_points);
    py_vertex_ids.resize({(int) number_of_corner_points});
    int* py_vertex_ids_ptr = static_cast<int*>(py_vertex_ids.request().ptr);
    for (std::size_t i_ctps{}; i_ctps < number_of_corner_points; i_ctps++) {
      py_vertex_ids_ptr[i_ctps] = static_cast<int>(vertex_ids[i_ctps]);
    }

    // Edges
    assert(edge_information.size() > 0);
    py::array_t<int> py_edges = py::array_t<int>(edge_information.size() * 3);
    py_edges.resize({(int) edge_information.size(), (int) 3});
    int* py_edges_ptr = static_cast<int*>(py_edges.request().ptr);
    for (std::size_t i_edge{}; i_edge < edge_information.size(); i_edge++) {
      for (std::size_t i_face{}; i_face < 3ul; i_face++) {
        py_edges_ptr[i_edge * 3ul + i_face] =
            static_cast<int>(edge_information[i_edge][i_face]);
      }
    }

    // Boundaries
    assert(boundaries.size() > 0);
    py::array_t<int> py_boundaries = py::array_t<int>(boundaries.size() * 4);
    py_boundaries.resize({(int) boundaries.size(), (int) 4});
    int* py_boundaries_ptr = static_cast<int*>(py_boundaries.request().ptr);
    for (std::size_t i_boundary{}; i_boundary < boundaries.size();
         i_boundary++) {
      for (std::size_t i_id{}; i_id < 4ul; i_id++) {
        py_boundaries_ptr[i_boundary * 4ul + i_id] =
            static_cast<int>(boundaries[i_boundary][i_id]);
      }
    }

    return py::make_tuple(py_connectivity,
                          py_vertex_ids,
                          py_edges,
                          py_boundaries,
                          is_structured);
  } else {
    throw std::runtime_error("Dimension mismatch");
  }
}

/**
 * @brief Extract all Boundary Patches and store them in a python list
 *
 * @param spline_list List of splines
 * @param interfaces interfaces, with negative values for boundary elements
 * @return py::list
 */
py::list ExtractAllBoundarySplines(const py::list& spline_list,
                                   const py::array_t<int>& interfaces,
                                   const int& n_threads) {
  // Check input data
  if (static_cast<int>(py::len(spline_list)) != interfaces.shape(0)) {
    splinepy::utils::PrintAndThrowError(
        "Number of splines in list (",
        py::len(spline_list),
        ") and number of elements in interfaces (",
        interfaces.shape(0),
        ") does not match.");
  }

  if (n_threads < 1) {
    splinepy::utils::PrintAndThrowError(
        "Number of threads must be positive integer.");
  }
  // Auxiliary data
  py::list boundary_splines{};
  std::vector<py::list> lists_to_concatenate(n_threads);
  const int* interface_ptr = static_cast<int*>(interfaces.request().ptr);
  const int n_patches = interfaces.shape(0);
  const int n_faces = interfaces.shape(1);
  const int para_dim_ = n_faces / 2;
  const auto cpp_spline_list =
      ListOfPySplinesToVectorOfCoreSplines(spline_list);
  const int chunk_size = std::div((n_patches + n_threads - 1), n_threads).quot;

  // This approach is a work-around for parallel execution
  auto extract_boundaries = [&](const int start, const int) {
    // start : process-ID
    // end : unused hence no referencing

    // Auxiliary variables
    const int start_index = start * chunk_size;
    const int end_index = (start + 1) * chunk_size;
    auto& boundaries_local = lists_to_concatenate[start];
    // Start extraction (remaining order)
    for (int i{start_index}; i < end_index; i++) {
      for (int j{}; j < n_faces; j++) {
        if (interface_ptr[i * n_faces + j] < 0) {
          boundaries_local.append(
              // Extract Boundary splines
              PySpline(cpp_spline_list[i]->SplinepyExtractBoundary(j)));
        }
      }
    }
  };

  // Execute in parallel
  splinepy::utils::NThreadExecution(extract_boundaries, n_threads, n_threads);

  // Concatenate list of boundaries - should only copy pointers
  for (auto& entries : lists_to_concatenate) {
    boundary_splines += entries;
  }

  return boundary_splines;
}

/**
 * @brief  Adds a Boundary using a seed and G-continuity on boundary-splines
 *
 * This function might be a slight overkill, as it assignes all functions an ID,
 * even when previously assigned a different ID -> Future Project
 *
 * @param boundary_splines boundary patches
 * @param boundary_interfaces interfaces between boundary splines
 * @param global_interfaces global interfaces (in between "volume"-patches)
 * @param tolerance tolerance to be considered g1 (1 - cos(phi) < tolerance)
 * @param n_threads number of threads for parallel processing
 * @return int number of new boundaries
 */
int AddBoundariesFromContinuity(const py::list& boundary_splines,
                                const py::array_t<int>& boundary_interfaces,
                                py::array_t<int>& global_interfaces,
                                const double& tolerance,
                                const int& n_threads) {
  // Check input data
  if (static_cast<int>(py::len(boundary_splines))
      != boundary_interfaces.shape(0)) {
    splinepy::utils::PrintAndThrowError(
        "Number of splines in list (",
        py::len(boundary_splines),
        ") and number of elements in connectivity (",
        boundary_interfaces.shape(0),
        ") does not match.");
  }

  // Provide auxiliary values
  const auto cpp_spline_list =
      ListOfPySplinesToVectorOfCoreSplines(boundary_splines);
  const int n_boundary_patches{static_cast<int>(boundary_interfaces.shape(0))};
  const int n_faces_per_boundary_patch{
      static_cast<int>(boundary_interfaces.shape(1))};
  const int para_dim_{n_faces_per_boundary_patch / 2};
  const int dim_ = cpp_spline_list[0]->SplinepyDim();
  const int* boundary_interfaces_ptr =
      static_cast<int*>(boundary_interfaces.request().ptr);
  int* global_interfaces_ptr =
      static_cast<int*>(global_interfaces.request().ptr);

  // Auxiliary Lambdas to keep code clean
  // Check if to tangential vectors are g1 (tol > cos(phi))
  auto areG1 = [&tolerance, &dim_](const std::vector<double>& vec0,
                                   const std::vector<double>& vec1) -> bool {
    // Checks in Debug
    assert(static_cast<int>(vec0.size()) == dim_);
    assert(static_cast<int>(vec1.size()) == dim_);

    // Start actual computation
    double norm0{}, norm1{}, dot_p{};
    for (int i{}; i < dim_; i++) {
      norm0 += vec0[i] * vec0[i];
      norm1 += vec1[i] * vec1[i];
      dot_p += vec0[i] * vec1[i];
    }
    return (tolerance > abs(1 - abs(dot_p) / std::sqrt(norm0 * norm1)));
  };

  // Identify face_id in adjacent patch
  auto face_id_in_neighbor =
      [&boundary_interfaces_ptr,
       &n_faces_per_boundary_patch](const int& base_patch_id,
                                    const int& neighbor_patch_id) -> int {
    // Loop over adjacent elements until base_patch_id is found
    for (int i{}; i < n_faces_per_boundary_patch; i++) {
      if (boundary_interfaces_ptr[neighbor_patch_id * n_faces_per_boundary_patch
                                  + i]
          == base_patch_id) {
        return i;
      }
    }
    // This part should never be reached
    splinepy::utils::PrintAndThrowError("Interface connectivity has errors, "
                                        "unidirectional interface detected.");
    return -1;
  };

  // Tangential Vector on boundary based on its derivative
  auto tangential_vector = [&cpp_spline_list, &para_dim_, &dim_](
                               const int& patch_id,
                               const int& face_id) -> std::vector<double> {
    // init return value (are default initialized to 0)
    std::vector<double> para_coord(para_dim_), bounds(2 * para_dim_),
        tangential_vector(dim_);
    std::vector<int> orders(para_dim_);

    // Auxiliary values
    const int axis_dim = face_id / 2;
    const int is_in_front = face_id % 2;

    // Parametric Bounds
    const auto& spline = cpp_spline_list[patch_id];
    spline->SplinepyParametricBounds(bounds.data());

    for (int i{}; i < para_dim_; i++) {
      if (i == axis_dim) {
        para_coord[i] = bounds[i + is_in_front * para_dim_];
        orders[i] = 1;
      } else {
        para_coord[i] = .5 * (bounds[i + para_dim_] + bounds[i]);
      }
    }
    spline->SplinepyDerivative(para_coord.data(),
                               orders.data(),
                               tangential_vector.data());
    return tangential_vector;
  };

  // Start Computations ------------------------------------------------ //
  // while the actual propagation needs to be performed in serial, the
  // precomputation of interface tolerances can be performed in parallel
  std::vector<bool> faces_are_g1(n_faces_per_boundary_patch
                                 * n_boundary_patches);
  auto precompute_tolerances = [&](const int start, const int end) {
    // Loop over relevant faces
    for (int i{start}; i < end; i++) {
      // Loop over faces
      for (int j{}; j < n_faces_per_boundary_patch; j++) {
        const int& adjacent_id =
            boundary_interfaces_ptr[i * n_faces_per_boundary_patch + j];

        if (adjacent_id < i) {
          // only compute if the adjacent neighbor has higher id to prevent
          // double the work
          continue;
        }
        // Get tangential vector of current patch
        const std::vector<double> vec0 = tangential_vector(i, j);

        // Get corresponding tangential vector of neighbor patch
        const int adjacent_face_id = face_id_in_neighbor(i, adjacent_id);
        const std::vector<double> vec1 =
            tangential_vector(adjacent_id, adjacent_face_id);

        // Check tolerance
        const bool is_g1 = areG1(vec0, vec1);
        faces_are_g1[i * n_faces_per_boundary_patch + j] = is_g1;
        faces_are_g1[adjacent_id * n_faces_per_boundary_patch
                     + adjacent_face_id] = is_g1;
      }
    }
  };

  // Execute in parallel
  splinepy::utils::NThreadExecution(precompute_tolerances,
                                    n_boundary_patches,
                                    n_threads);

  // std::vector can use less memory for bools
  std::vector<bool> is_assigned(n_boundary_patches); // defaults false
  std::vector<int> new_boundary_id(n_boundary_patches);
  std::vector<int> queued_splines{};

  // Start Assignement
  // Loop over all patches
  int current_max_id{1};
  for (int i{}; i < n_boundary_patches; i++) {
    if (is_assigned[i]) {
      continue;
    }
    new_boundary_id[i] = current_max_id;
    is_assigned[i] = true;
    queued_splines.push_back(i);

    // Start propagation
    while (!queued_splines.empty()) {
      const int current_id = queued_splines.back();
      queued_splines.pop_back();
      for (int i_face{}; i_face < n_faces_per_boundary_patch; i_face++) {
        const int combined_index =
            current_id * n_faces_per_boundary_patch + i_face;
        // Is the neighborface G1
        if (faces_are_g1[combined_index]) {
          const int& adjacent_id = boundary_interfaces_ptr[combined_index];
          // Check if the adjacent patch is already assigned
          if (is_assigned[adjacent_id]) {
            continue;
          } else {
            // Assign a BID and continue
            new_boundary_id[adjacent_id] = current_max_id;
            is_assigned[adjacent_id] = true;
            queued_splines.push_back(adjacent_id);
          }
        }
      }
    }

    // End propagation and increase id
    current_max_id++;
  }

  // Assign the new boundary ids to the old interface-vector
  const int& n_interfaces = global_interfaces.size();
  int counter{};
  for (int i{}; i < n_interfaces; i++) {
    if (global_interfaces_ptr[i] < 0) {
      global_interfaces_ptr[i] = -new_boundary_id[counter];
      counter++;
    }
  }
  if (counter != n_boundary_patches) {
    splinepy::utils::PrintAndThrowError(
        counter,
        " new boundary ids were assigned, however ",
        n_boundary_patches,
        " were expected, which means information was lost. Abort mission");
  }

  return current_max_id;
}

} // namespace splinepy::py
