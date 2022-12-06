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
#include <splinepy/utils/print.hpp>

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
SplineLibIoSpline PySplineToSplineLibIoSpline(PySpline& pyspline) {
  return std::dynamic_pointer_cast<typename SplineLibIoSpline::element_type>(
      splinepy::py::SameSplineWithKnotVectors(pyspline).Core());
}

/// convert list of PySplines to vector of splinelib SplineItems
SplineLibIoSplines ListOfPySplinesToSplineLibIoSplines(py::list pysplines) {
  // prepare return obj
  SplineLibIoSplines io_splines;
  io_splines.reserve(pysplines.size());
  for (py::handle pys : pysplines) {
    io_splines.emplace_back(
        PySplineToSplineLibIoSpline(py::cast<PySpline&>(pys)));
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
  using ResolutionsPerSpline = splinelib::sources::input_output::vtk::
      NumbersOfParametricCoordinates;
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
  const std::size_t problem_dimension = corner_vertex_buffer.shape[1];
  const std::size_t number_of_corner_points = corner_vertex_buffer.shape[0];
  // Check if mesh can be used for mfem mesh
  bool is_structured{true};
  if (problem_dimension == 2) {
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
      for (std::size_t i_dim{}; i_dim < problem_dimension; i_dim++) {
        vertex[i_dim] = corner_ptr[i_c * problem_dimension + i_dim];
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
    const auto [connectivity,
                vertex_ids,
                edge_information,
                boundaries] = [&]() {
      try {
        return bezman::utils::algorithms::ExtractMFEMInformation(
            corner_vertices,
            maxvertex - minvertex,
            tolerance);
      } catch (...) {
        is_structured = false;
        return std::make_tuple(
            // Connectivity
            bezman::utils::algorithms::FindConnectivity(corner_vertices,
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

  } else if (problem_dimension == 3) {
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
      for (std::size_t i_dim{}; i_dim < problem_dimension; i_dim++) {
        vertex[i_dim] = corner_ptr[i_c * problem_dimension + i_dim];
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
    const auto [connectivity,
                vertex_ids,
                edge_information,
                boundaries] = [&]() {
      try {
        return bezman::utils::algorithms::ExtractMFEMInformation(
            corner_vertices,
            maxvertex - minvertex,
            tolerance);
      } catch (...) {
        is_structured = false;
        return std::make_tuple(
            // Connectivity
            bezman::utils::algorithms::FindConnectivity(corner_vertices,
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

} // namespace splinepy::py
