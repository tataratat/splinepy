#pragma once

#include <cstdlib>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_map>

// pybind
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

//
#include "splinepy/py/py_spline.hpp"
#include "splinepy/splines/null_spline.hpp"
#include "splinepy/utils/print.hpp"

namespace splinepy::py {

namespace py = pybind11;

// alias for frequently used vectors
using CoreSplineVector =
    std::vector<std::shared_ptr<splinepy::splines::SplinepyBase>>;
using IntVector = splinepy::utils::DefaultInitializationVector<int>;
using IntVectorVector = splinepy::utils::DefaultInitializationVector<IntVector>;
using DoubleVector = splinepy::utils::DefaultInitializationVector<double>;

/// @brief Extracts CoreSpline from list of PySplines.
/// @param pysplines
/// @param nthreads
/// @return
inline std::vector<PySpline::CoreSpline_>
ToCoreSplineVector(py::list pysplines, const int nthreads = 1) {
  // prepare return obj
  const int n_splines = static_cast<int>(pysplines.size());
  std::vector<PySpline::CoreSpline_> core_splines(n_splines);

  auto to_core = [&](int begin, int end) {
    for (int i{begin}; i < end; ++i) {
      core_splines[i] =
          pysplines[i].template cast<std::shared_ptr<PySpline>>()->Core();
    }
  };
  splinepy::utils::NThreadExecution(to_core, n_splines, nthreads);

  return core_splines;
}

/// @brief Overload to allow None entries. Will be filled with null spline
/// @param pysplines
/// @param para_dim_if_none
/// @param dim_if_none
/// @param nthreads
/// @return
inline std::vector<PySpline::CoreSpline_>
ToCoreSplineVector(py::list pysplines,
                   const int para_dim_if_none,
                   const int dim_if_none,
                   const int nthreads) {
  // prepare return obj
  const int n_splines = static_cast<int>(pysplines.size());
  std::vector<PySpline::CoreSpline_> core_splines(n_splines);

  auto to_core = [&](int begin, int end) {
    for (int i{begin}; i < end; ++i) {
      // get accessor
      auto spline = pysplines[i];

      // if None, fill null spline
      if (spline.is_none()) {
        core_splines[i] =
            splinepy::splines::kNullSplineLookup[para_dim_if_none - 1]
                                                [dim_if_none - 1];
        continue;
      }
      core_splines[i] =
          spline.template cast<std::shared_ptr<PySpline>>()->Core();
    }
  };
  splinepy::utils::NThreadExecution(to_core, n_splines, nthreads);

  return core_splines;
}

/// @brief Converts a CoreSplineVector to a Python spline list
/// @param splist
/// @return py::list
inline py::list ToPySplineList(CoreSplineVector& splist) {
  // prepare return obj
  const int n_splines = static_cast<int>(splist.size());

  // to return
  py::list pyspline_list(n_splines);

  // cast - single threads, as we should respect GIL
  for (int i{}; i < n_splines; ++i) {
    pyspline_list[i] = std::make_shared<PySpline>(splist[i])->ToDerived();
  }

  return pyspline_list;
}

/*
 * Sort Vector using lambda expressions
 * https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
 */
template<typename T>
std::vector<std::size_t> IndexListSort(const std::vector<T>& v) {
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {
    return v[i1] < v[i2];
  });
  return idx;
}

/// @brief raises if there's any mismatch between specified properties and all
/// the entries in the vector.
/// @param splist
/// @param name
/// @param para_dim
/// @param dim
/// @param degrees
/// @param control_mesh_resolutions
/// @param nthreads
inline void RaiseMismatch(const CoreSplineVector& splist,
                          const std::string name,
                          const int para_dim,
                          const int dim,
                          const IntVector& degrees,
                          const IntVector& control_mesh_resolutions,
                          const int nthreads) {
  // for verbose output
  std::unordered_map<std::string, IntVectorVector> mismatches{};

  // check flags
  bool check_name{}, check_para_dim{}, check_dim{}, check_degrees{},
      check_control_mesh_resolutions{};
  // input vector sizes
  const int d_size = degrees.size();
  const int cmr_size = control_mesh_resolutions.size();
  // true para_dim to use as reference -> will be updated / compared
  int ref_para_dim{para_dim};

  // parse input and allocate vector for mismatch book keeping
  if (name.size() > 0) {
    check_name = true;
    mismatches["name"].resize(nthreads);
  }
  if (dim > 0) {
    check_dim = true;
    mismatches["dim"].resize(nthreads);
  }
  if (d_size > 0) {
    check_degrees = true;
    ref_para_dim = d_size;
    mismatches["degrees"].resize(nthreads);
  }
  if (cmr_size > 0) {
    check_control_mesh_resolutions = true;
    mismatches["control_mesh_resolutions"].resize(nthreads);
  }
  if (para_dim > 0 || check_degrees || check_control_mesh_resolutions) {
    check_para_dim = true;

    // sanity check
    for (const auto& candidate : {para_dim, d_size, cmr_size}) {
      if (candidate != ref_para_dim && candidate > 0) {
        splinepy::utils::PrintAndThrowError(
            "Mismatch in given para_dim (",
            para_dim,
            "), size of degrees, (",
            d_size,
            "), size of control_mesh_resolutions (",
            cmr_size,
            ").");
      }
    }
    mismatches["para_dim"].resize(nthreads);
  }

  // lambda for nthread comparison
  auto check_mismatch_step = [&](int begin, int total_) {
    // in step-style query, begin is i_thread
    const int thread_index = begin;
    // alloc vectors incase we need to compare
    IntVector spline_degree(ref_para_dim), spline_cmr(ref_para_dim);

    for (int i{begin}; i < total_; i += nthreads) {
      // get spline to check
      const auto& spline = splist[i];

      // skip null splines
      if (spline->SplinepyIsNull()) {
        continue;
      }

      // default value for para_dim match
      bool para_dim_matches{true};

      // check
      if (check_name && (name != spline->SplinepySplineName())) {
        mismatches["name"][thread_index].push_back(i);
      }
      if (check_para_dim && (ref_para_dim != spline->SplinepyParaDim())) {
        mismatches["para_dim"][thread_index].push_back(i);
        para_dim_matches = false;
      }
      if (check_dim && (dim != spline->SplinepyDim())) {
        mismatches["dim"][thread_index].push_back(i);
      }
      // check properties that are relevent iff para_dim matches
      if (para_dim_matches) {
        if (check_degrees) {
          spline->SplinepyCurrentProperties(spline_degree.data(),
                                            nullptr,
                                            nullptr,
                                            nullptr);
          if (spline_degree != degrees) {
            mismatches["degrees"][thread_index].push_back(i);
          }
        }
        if (check_control_mesh_resolutions) {
          spline->SplinepyControlMeshResolutions(spline_cmr.data());
          if (spline_cmr != control_mesh_resolutions) {
            mismatches["control_mesh_resolutions"][thread_index].push_back(i);
          }
        }
      }
    }
  };

  splinepy::utils::NThreadExecution(check_mismatch_step,
                                    static_cast<int>(splist.size()),
                                    nthreads,
                                    splinepy::utils::NThreadQueryType::Step);

  // prepare output or maybe exit
  std::unordered_map<std::string, IntVector> concat_mismatches{};
  bool raise{false};

  for (const auto& [key, mismatch_per_threads] : mismatches) {
    auto& concat_vec = concat_mismatches[key];
    for (const auto& mismatch : mismatch_per_threads) {
      const auto m_size = mismatch.size();
      if (m_size != 0) {
        concat_vec.insert(concat_vec.end(), mismatch.begin(), mismatch.end());
      }
    }
    if (concat_vec.size() != 0) {
      raise = true;
    }
  };

  // everything matches
  if (!raise) {
    return;
  }

  // form mismatch info
  std::string mismatch_info{};
  for (const auto& [key, concat_mismatch] : concat_mismatches) {
    mismatch_info += "\n[" + key + "] : ";
    for (const auto& ids : concat_mismatch) {
      mismatch_info += std::to_string(ids) + ", ";
    }
  }

  splinepy::utils::PrintAndThrowError("Found mismatches.", mismatch_info);
}

/// @brief raises if there's any mismatch between two given vectors of splines.
/// @param splist0
/// @param splist1
/// @param name
/// @param para_dim
/// @param dim
/// @param degrees
/// @param control_mesh_resolutions
/// @param nthreads
inline void RaiseMismatch(const CoreSplineVector& splist0,
                          const CoreSplineVector& splist1,
                          const bool name,
                          const bool para_dim,
                          const bool dim,
                          const bool degrees,
                          const bool control_mesh_resolutions,
                          const int nthreads) {
  // len check
  if (splist0.size() != splist1.size()) {
    splinepy::utils::PrintAndThrowError(
        "Size mismatch between first and second list of splines.",
        splist0.size(),
        "vs",
        splist1.size());
  }
  // for verbose output
  std::unordered_map<std::string, IntVectorVector> mismatches{};
  // parse input and allocate vector for mismatch book keeping
  if (name) {
    mismatches["name"].resize(nthreads);
  }
  if (dim) {
    mismatches["dim"].resize(nthreads);
  }
  if (degrees) {
    mismatches["degrees"].resize(nthreads);
  }
  if (control_mesh_resolutions) {
    mismatches["control_mesh_resolutions"].resize(nthreads);
  }
  if (para_dim || degrees || control_mesh_resolutions) {
    mismatches["para_dim"].resize(nthreads);
  }

  // lambda for nthread comparison
  auto check_mismatch_step = [&](int begin, int total_) {
    // in step-style query, begin is i_thread
    const int thread_index = begin;

    // alloc some tmp vector if needed
    IntVector int_vec0, int_vec1;

    for (int i{begin}; i < total_; i += nthreads) {
      // get spline to check
      const auto& spline0 = splist0[i];
      const auto& spline1 = splist1[i];

      // default value for para_dim match
      bool para_dim_matches{true};

      // check dims
      if (para_dim
          && (spline0->SplinepyParaDim() != spline1->SplinepyParaDim())) {
        mismatches["para_dim"][thread_index].push_back(i);
        para_dim_matches = false;
      }
      if (dim && (spline0->SplinepyDim() != spline1->SplinepyDim())) {
        mismatches["dim"][thread_index].push_back(i);
      }
      // null spline should have at least same dims. but nothing else matters.
      // skip null splines
      if (spline1->SplinepyIsNull() || spline0->SplinepyIsNull()) {
        continue;
      }

      if (name
          && (spline0->SplinepySplineName() != spline1->SplinepySplineName())) {
        mismatches["name"][thread_index].push_back(i);
      }

      // check properties that are relevent iff para_dim matches
      if (para_dim_matches && (degrees || control_mesh_resolutions)) {
        // alloc some space
        int_vec0.resize(spline0->SplinepyParaDim());
        int_vec1.resize(spline1->SplinepyParaDim());

        if (degrees) {
          spline0->SplinepyCurrentProperties(int_vec0.data(),
                                             nullptr,
                                             nullptr,
                                             nullptr);
          spline1->SplinepyCurrentProperties(int_vec1.data(),
                                             nullptr,
                                             nullptr,
                                             nullptr);

          if (int_vec0 != int_vec1) {
            mismatches["degrees"][thread_index].push_back(i);
          }
        }

        if (control_mesh_resolutions) {
          spline0->SplinepyControlMeshResolutions(int_vec0.data());
          spline1->SplinepyControlMeshResolutions(int_vec1.data());

          if (int_vec0 != int_vec1) {
            mismatches["control_mesh_resolutions"][thread_index].push_back(i);
          }
        }
      }
    }
  };

  splinepy::utils::NThreadExecution(check_mismatch_step,
                                    static_cast<int>(splist0.size()),
                                    nthreads,
                                    splinepy::utils::NThreadQueryType::Step);

  // prepare output or maybe exit
  std::unordered_map<std::string, IntVector> concat_mismatches{};
  bool raise{false};

  for (const auto& [key, mismatch_per_threads] : mismatches) {
    auto& concat_vec = concat_mismatches[key];
    for (const auto& mismatch : mismatch_per_threads) {
      const auto m_size = mismatch.size();
      if (m_size != 0) {
        concat_vec.insert(concat_vec.end(), mismatch.begin(), mismatch.end());
      }
    }
    if (concat_vec.size() != 0) {
      raise = true;
    }
  };

  // everything matches
  if (!raise) {
    return;
  }

  // form mismatch info
  std::string mismatch_info{};
  for (const auto& [key, concat_mismatch] : concat_mismatches) {
    mismatch_info += "\n[" + key + "] : ";
    for (const auto& ids : concat_mismatch) {
      mismatch_info += std::to_string(ids) + ", ";
    }
  }

  splinepy::utils::PrintAndThrowError("Found mismatches.", mismatch_info);
}

/**
 * @brief Determine the connectivity from center-vertices, assuming nothing of
 * the underlying grid
 *
 * Duplicate Points are not eliminated, assuming that a maximum of two points
 * are equivalent. If this is not the case an exception is thrown. In theory
 * this has complexity O(nlogn) whereas a KDTree has complexity O(n (logn)^dim).
 *
 * @tparam PhysicalPointType Type of Point coordinates
 * @tparam ScalarType Type determining the precision
 * @tparam parametric_dimension dimension of the object (e.g. surface in 3D)
 * @tparam boolean check_orientation to check if neighboring elements match
 *                          structured grid
 * @param face_center_vertices vertices in the centers of spline-surfaces
 * @param metric used for preordering the vertices along a line
 * @param tolerance tolerance (distance between two vertices that are joined)

 * @return connectivity as a std::vector<std::array<...>>
 */
py::array_t<int>
FindConnectivityFromCenters(const py::array_t<double>& face_center_vertices,
                            const int parametric_dimension,
                            const py::array_t<double>& metric,
                            const double tolerance) {

  // -- Auxiliary data --
  const int number_of_element_faces = parametric_dimension * 2;
  const int number_of_patches =
      face_center_vertices.shape(0) / number_of_element_faces;
  const int physical_dimension = face_center_vertices.shape(1);
  const int number_of_center_vertices{
      static_cast<int>(face_center_vertices.shape(0))};
  const double tolerance_squared{tolerance * tolerance};
  const double* metric_ptr = static_cast<double*>(metric.request().ptr);
  const double* face_center_vertices_ptr =
      static_cast<double*>(face_center_vertices.request().ptr);

  // Consistency check
  if (!(face_center_vertices.shape(0) % number_of_element_faces == 0)) {
    splinepy::utils::PrintAndThrowError(
        "Number of corner vertices invalid. Must be a multiple of the number "
        "of vertices per patch");
  }
  if (!(metric.size() == face_center_vertices.shape(1))) {
    splinepy::utils::PrintAndThrowError(
        "Incompatible size for metric. Must match physical dimension of face "
        "center vertices");
  }
  if (!(number_of_center_vertices % number_of_element_faces == 0)) {
    splinepy::utils::PrintAndThrowError(
        "Inconsistent number of Center vertices. Must be divisible by "
        "parametric_dimension*2");
  }

  // Init return type
  py::array_t<int> connectivity(number_of_patches * number_of_element_faces);
  int* connectivity_ptr = static_cast<int*>(connectivity.request().ptr);
  // Init connectivity and metric value
  // (-1 : boundary, -2 : untouched)
  std::fill(connectivity_ptr,
            connectivity_ptr + number_of_patches * number_of_element_faces,
            -2);

  // Assure Metric is normed and non-zero
  const std::vector<double>& normed_metric = [&metric_ptr,
                                              &physical_dimension]() {
    double euclidian_norm = 0.;
    for (int i_phys{}; i_phys < physical_dimension; i_phys++) {
      euclidian_norm += metric_ptr[i_phys] * metric_ptr[i_phys];
    }
    std::vector<double> metric(physical_dimension, 1.);
    euclidian_norm = std::sqrt(euclidian_norm);
    if (euclidian_norm < 1e-20) {
      return metric;
    } else {
      const double inv_euclidian_norm = 1 / euclidian_norm;
      for (int i_phys{}; i_phys < physical_dimension; i_phys++) {
        metric[i_phys] = metric_ptr[i_phys] * inv_euclidian_norm;
      }
    }
    return metric;
  }();

  // Auxiliary function to determine distance between two points in face_center
  auto squared_euclidian_distance =
      [&face_center_vertices_ptr,
       &physical_dimension](const int& i_start, const int& i_end) -> double {
    double squared_euclidian_distance_{};
    for (int i_phys{}; i_phys < physical_dimension; i_phys++) {
      const double distance_c =
          face_center_vertices_ptr[i_end * physical_dimension + i_phys]
          - face_center_vertices_ptr[i_start * physical_dimension + i_phys];
      squared_euclidian_distance_ += distance_c * distance_c;
    }
    return squared_euclidian_distance_;
  };

  std::vector<double> scalar_metric{};
  scalar_metric.reserve(number_of_center_vertices);

  // Check Metric Dimension and Vector Size
  for (int i_vertex{}; i_vertex < number_of_center_vertices; i_vertex++) {
    double metric_v = normed_metric[0]
                      * face_center_vertices_ptr[i_vertex * physical_dimension];
    for (int j_phys{1}; j_phys < physical_dimension; j_phys++) {
      metric_v =
          normed_metric[j_phys]
          * face_center_vertices_ptr[i_vertex * physical_dimension + j_phys];
    }
    scalar_metric.push_back(metric_v);
  }

  // Sort Metric Vector
  const auto metric_order_indices = IndexListSort(scalar_metric);

  // Loop over points
  for (int lower_limit = 0; lower_limit < number_of_center_vertices - 1;
       lower_limit++) {
    // Loop over all points regardless of whether they have been touched or not,
    // and then check the validity of the connection Point already processed
    bool found_duplicate = false;
    // Now check allowed range for duplicates
    int upper_limit = lower_limit + 1;
    while (upper_limit < number_of_center_vertices
           && (scalar_metric[metric_order_indices[upper_limit]]
               - scalar_metric[metric_order_indices[lower_limit]])
                  < tolerance) {
      // Check if the two points are duplicates
      found_duplicate =
          squared_euclidian_distance(metric_order_indices[lower_limit],
                                     metric_order_indices[upper_limit])
          < tolerance_squared;
      if (found_duplicate) {
        break;
      } else {
        upper_limit++;
      }
    }

    // Now we have to check if the connection is valid
    // 1. If another connection is found, that means, that the point it connects
    //    to has a higher index in the metric tensor. If the current point does
    //    already have a neighbor, that means that more than one point connect
    //    -> Error
    if (found_duplicate) {
      // Check 1. (@todo EXCEPTION)
      if (connectivity_ptr[metric_order_indices[lower_limit]] != -2) {
        splinepy::utils::PrintAndThrowError(
            "Found conflicting interceptions, where more than two points are "
            "in the same position. Expected -2 got",
            connectivity_ptr[metric_order_indices[lower_limit]]);
      }

      // If both tests passed, update connectivity
      connectivity_ptr[metric_order_indices[lower_limit]] =
          static_cast<int>(metric_order_indices[upper_limit])
          / number_of_element_faces;
      connectivity_ptr[metric_order_indices[upper_limit]] =
          static_cast<int>(metric_order_indices[lower_limit])
          / number_of_element_faces;
    } else {
      // set Boundary-ID
      if (connectivity_ptr[metric_order_indices[lower_limit]] == -2) {
        connectivity_ptr[metric_order_indices[lower_limit]] = -1;
      }
    }
  }

  // Treat last remaining point in scalar metric vector
  if (connectivity_ptr[metric_order_indices[number_of_center_vertices - 1]]
      == -2) {
    connectivity_ptr[metric_order_indices[number_of_center_vertices - 1]] = -1;
  }

  // Resize buffer
  connectivity.resize({number_of_patches, number_of_element_faces});
  return connectivity;
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
  // Determine Metric
  const int physical_dimension = py_center_vertices.shape(1);
  const int number_of_corner_points = py_center_vertices.shape(0);
  py::array_t<double> maxvertex(physical_dimension),
      minvertex{physical_dimension};

  // Get pointers
  double* maxvertex_ptr = static_cast<double*>(maxvertex.request().ptr);
  double* minvertex_ptr = static_cast<double*>(minvertex.request().ptr);
  const double* py_center_vertices_ptr =
      static_cast<double*>(py_center_vertices.request().ptr);

  // Calculate limits of center vertices
  for (int i_dim{}; i_dim < physical_dimension; i_dim++) {
    minvertex_ptr[i_dim] = py_center_vertices_ptr[i_dim];
    maxvertex_ptr[i_dim] = py_center_vertices_ptr[i_dim];
  }
  for (int i_c{1}; i_c < number_of_corner_points; i_c++) {
    for (int i_dim{}; i_dim < physical_dimension; i_dim++) {
      maxvertex_ptr[i_dim] =
          std::max(py_center_vertices_ptr[i_c * physical_dimension + i_dim],
                   maxvertex_ptr[i_dim]);
      minvertex_ptr[i_dim] =
          std::min(py_center_vertices_ptr[i_c * physical_dimension + i_dim],
                   minvertex_ptr[i_dim]);
    }
  }

  // Determine Interfaces
  return FindConnectivityFromCenters(py_center_vertices,
                                     parametric_dimension,
                                     maxvertex - minvertex,
                                     tolerance);
}

/// @brief Orientation between two adjacent splines
///
/// If two splines share the same boundary this function retrieves their
/// orientation, by mapping the mappings of the parametric axis onto each other.
/// This is (among others) required for Gismo and Nutils export
///
/// @param pyspline_start Spline object from start
/// @param boundary_start Boundary ID from start spline
/// @param pyspline_end Spline object from end///to which is mapped
/// @param boundary_end Boundary ID of adjacent spline
/// @param tolerance Tolerance
/// @param int_mappings_ptr (output) integer mappings
/// @param bool_orientations_ptr (output) axis alignement
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

/// @brief Get the Boundary Orientations object
///
/// @param spline_list
/// @param base_id
/// @param base_face_id
/// @param neighbor_id
/// @param neighbor_face_id
/// @param tolerance
/// @param n_threads
/// @return py::tuple
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
  const auto cpp_spline_list = ToCoreSplineVector(spline_list);
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

/// @brief Extract all Boundary Patches and store them in a python list
///
/// @param spline_list List of splines
/// @param interfaces interfaces, with negative values for boundary elements
/// @param n_threads Number of threads
/// @return py::list
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
  const auto cpp_spline_list = ToCoreSplineVector(spline_list);
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
  const auto cpp_spline_list = ToCoreSplineVector(boundary_splines);

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

class PyMultiPatch;

/// @brief ToDerived for PyMultiPatches
/// @param core_obj
/// @return
py::object ToDerived(std::shared_ptr<PyMultiPatch> core_obj) {
  const auto to_derived = py::module_::import("splinepy").attr("to_derived");
  return to_derived(py::cast(core_obj));
}

/// @brief Multi-patch splines. Here, use of the words
///        "patch" and "spline" are interchangeable
class PyMultiPatch {
public:
  using CorePatches_ = std::vector<typename PySpline::CoreSpline_>;

  /// individual patches, a.k.a. splines
  py::list patches_;

  /// casted, core splines (patches) from patches_ (list of PySpline).
  CorePatches_ core_patches_;

  /// sub patches - similar concept to subelement / face elements
  std::shared_ptr<PyMultiPatch> sub_multi_patch_ = nullptr;

  /// center of each sub patch.
  py::array_t<double> sub_patch_centers_;

  /// ids for boundary_patches from sub_patches_
  py::array_t<int> boundary_patch_ids_;

  /// Boundary multi patches. Should consist of one smaller parametric dim
  std::shared_ptr<PyMultiPatch> boundary_multi_patch_ = nullptr;

  /// patch-to-patch connectivity information
  /// shape: (n_patches, n_boundary_elements)
  py::array_t<int> interfaces_;

  /// shape: (n_patches,)
  py::array_t<int> boundary_ids_;

  /// @brief Fields - raw list format
  py::list field_list_;

  /// @brief Fields - they are saved as multi-patches
  std::vector<std::shared_ptr<PyMultiPatch>> field_multi_patches_;

  /// default number of threads for all the operations besides queries
  int n_default_threads_{1};

  /// hint flag that can be used to accelerate property decisions
  bool same_parametric_bounds_{false};

  /// internal flag to indicate that this a field multipatch that may contain
  /// null spline
  bool has_null_splines_{false};

  /// default tolerance
  double tolerance_{1e-11};

  /// ctor
  PyMultiPatch() = default;

  /// move ctor
  PyMultiPatch(PyMultiPatch&& other) noexcept = default;

  /// @brief list (iterable) init -> pybind will cast to list for us if needed
  /// @param patches
  /// @param n_default_threads
  /// @param same_parametric_bounds
  PyMultiPatch(py::list& patches,
               const int n_default_threads = 1,
               const bool same_parametric_bounds = false)
      : n_default_threads_(n_default_threads),
        same_parametric_bounds_(same_parametric_bounds) {
    SetPatchesNThreads(patches, n_default_threads);
  }

  /// @brief Internal use only. used to save field splines.
  /// @param patches
  /// @param n_default_threads
  /// @param para_dim_if_none
  /// @param dim_if_none
  PyMultiPatch(py::list& patches,
               const int n_default_threads,
               const int para_dim_if_none,
               const int dim_if_none)
      : n_default_threads_(n_default_threads),
        has_null_splines_(true) {
    SetPatchesWithNullSplines(patches,
                              n_default_threads,
                              para_dim_if_none,
                              dim_if_none);
  }

  /// @brief given other shared pointer, tries to steal the contents with move
  /// ctor. Expected use is only for  to_derived call.
  /// @param other
  PyMultiPatch(std::shared_ptr<PyMultiPatch>& other)
      : PyMultiPatch(std::move(*other)) {}

  /// @brief Gets core patches
  CorePatches_& CorePatches() {
    if (core_patches_.size() == 0) {
      splinepy::utils::PrintAndThrowError("No splines/patches set");
    }
    return core_patches_;
  }

  /// @brief Clears
  void Clear() {
    patches_ = py::list();
    core_patches_ = CorePatches_{};
    sub_multi_patch_ = nullptr;
    boundary_patch_ids_ = py::array_t<int>();
    boundary_multi_patch_ = nullptr;
    interfaces_ = py::array_t<int>();
    boundary_ids_ = py::array_t<int>();
    sub_patch_centers_ = py::array_t<double>();
    field_list_ = py::list();
    field_multi_patches_ = {};
  }

  /// @brief returns para_dim of the first patch
  /// @return
  int ParaDim() { return CorePatches()[0]->SplinepyParaDim(); }

  /// @brief returns dim of the first patch
  /// @return
  int Dim() { return CorePatches()[0]->SplinepyDim(); }

  /// @brief my name is MultiPatch
  /// @return
  std::string Name() { return "MultiPatch"; }

  /// @brief MultiPatch with n patches
  /// @return
  std::string WhatAmI() {
    return "MultiPatch with " + std::to_string(core_patches_.size())
           + " patches.";
  }

  /// could overload, but it won't be clear for pybind
  /// @brief default patches setter, exposed to python
  /// @param patches
  void SetPatchesDefault(py::list& patches) {
    SetPatchesNThreads(patches, n_default_threads_);
  }

  /// @brief patches setter with nthreads
  /// @param patches
  /// @param nthreads
  void SetPatchesNThreads(py::list& patches, const int nthreads = 1) {
    // clear first, incase this is a new setting
    Clear();

    // register patches
    patches_ = patches;

    // convert
    core_patches_ = ToCoreSplineVector(patches, nthreads);

    // check only dim mismatch
    RaiseMismatch(core_patches_,
                  "",
                  core_patches_[0]->SplinepyParaDim(),
                  core_patches_[0]->SplinepyDim(),
                  {},
                  {},
                  nthreads);
  }

  /// @brief patches setter for field patches that may contain null splines.
  /// @param patches
  /// @param nthreads
  /// @param para_dim_if_none
  /// @param dim_if_none
  void SetPatchesWithNullSplines(py::list& patches,
                                 const int nthreads,
                                 const int para_dim_if_none,
                                 const int dim_if_none) {
    // clear
    Clear();

    // save
    patches_ = patches;

    // convert using null spline friendly converter
    core_patches_ =
        ToCoreSplineVector(patches, para_dim_if_none, dim_if_none, nthreads);

    // check only dim mismatch
    RaiseMismatch(core_patches_,
                  "",
                  core_patches_[0]->SplinepyParaDim(),
                  core_patches_[0]->SplinepyDim(),
                  {},
                  {},
                  nthreads);
  }
  /// @brief Gets patches
  py::list GetPatches() { return patches_; }

  /// @brief returns sub patches as multi patches, uses default nthreads
  /// @return
  std::shared_ptr<PyMultiPatch> SubMultiPatch() {
    // quick exit if they exist.
    if (sub_multi_patch_) {
      return sub_multi_patch_;
    }

    const int n_splines = CorePatches().size();
    // to accumulate
    const int n_boundary = ParaDim() * 2;
    const int n_boundaries = n_splines * n_boundary;

    // prepare output
    CoreSplineVector out_boundaries(n_boundaries);

    // prepare lambda
    auto boundary_extract = [&](int begin, int end) {
      for (int i{begin}; i < end; ++i) {
        // start of the offset
        const int offset = i * n_boundary;
        // deref this spline
        auto& spline = *core_patches_[i];
        for (int j{}; j < n_boundary; ++j) {
          out_boundaries[offset + j] = spline.SplinepyExtractBoundary(j);
        }
      }
    };

    splinepy::utils::NThreadExecution(boundary_extract,
                                      n_splines,
                                      n_default_threads_);

    // create multi patch to return
    sub_multi_patch_ = std::make_shared<PyMultiPatch>();

    // set both core and py splines
    sub_multi_patch_->patches_ = ToPySplineList(out_boundaries);
    sub_multi_patch_->core_patches_ = std::move(out_boundaries);

    return sub_multi_patch_;
  }

  /// @brief returns derived class of SubMultiPatch
  /// @return
  py::object PySubMultiPatch() { return ToDerived(SubMultiPatch()); }

  /// @brief Computes centers of subpatches
  /// @return
  py::array_t<double> SubPatchCenters() {
    // return saved
    if (sub_patch_centers_.size() > 0) {
      return sub_patch_centers_;
    }

    // prepare output
    // from here we assume that all the splines have the same para_dim and dim
    const int n_splines = CorePatches().size();
    const int para_dim = ParaDim();
    const int dim = Dim();
    const int n_queries = 2 * para_dim;
    const int n_total = n_queries * n_splines;
    sub_patch_centers_ = py::array_t<double>({n_total, dim});
    double* sub_patch_centers_ptr =
        static_cast<double*>(sub_patch_centers_.request().ptr);

    // pre-compute boundary centers
    DoubleVector para_bounds;
    double* para_bounds_ptr;
    if (!same_parametric_bounds_) {
      para_bounds.resize(n_total * para_dim);
      para_bounds_ptr = para_bounds.data();
      const int stride = n_queries * para_dim;

      auto calc_para_bounds = [&](int begin, int end) {
        for (int i{begin}; i < end; ++i) {
          splinepy::splines::helpers::ScalarTypeBoundaryCenters(
              *core_patches_[i],
              &para_bounds_ptr[stride * i]);
        }
      };

      // exe
      splinepy::utils::NThreadExecution(calc_para_bounds,
                                        n_splines,
                                        n_default_threads_);
    }

    auto calc_sub_patch_centers_step = [&](int begin, int total_) {
      // each thread needs one query
      DoubleVector queries_vector; /* unused if same_parametric_bounds=true*/
      double* queries;

      // pre compute boundary centers if para bounds are the same
      if (same_parametric_bounds_) {
        queries_vector.resize(2 * para_dim * para_dim);
        queries = queries_vector.data();
        splinepy::splines::helpers::ScalarTypeBoundaryCenters(*core_patches_[0],
                                                              queries);
      }

      for (int i{begin}; i < total_; i += n_default_threads_) {
        const auto [i_spline, i_query] = std::div(i, n_queries);

        // get ptr start
        if (!same_parametric_bounds_) {
          queries = &para_bounds_ptr[i_spline * n_queries * para_dim];
        }

        // eval
        core_patches_[i_spline]->SplinepyEvaluate(
            &queries[i_query * para_dim],
            &sub_patch_centers_ptr[(i_spline * n_queries + i_query) * dim]);
      }
    };

    splinepy::utils::NThreadExecution(calc_sub_patch_centers_step,
                                      n_total,
                                      n_default_threads_,
                                      splinepy::utils::NThreadQueryType::Step);

    return sub_patch_centers_;
  }

  /// @brief setter and getter for interface info. runs array shape check for
  /// setter
  /// @param interfaces
  py::array_t<int> Interfaces(const py::array_t<int>& interfaces) {
    // check if we should set or get
    // set
    if (interfaces.size() > 0) {
      // size check
      splinepy::py::CheckPyArrayShape(
          interfaces,
          {static_cast<int>(CorePatches().size()),
           static_cast<int>(CorePatches()[0]->SplinepyParaDim() * 2)},
          true);

      interfaces_ = interfaces;
    } else if (interfaces_.size() == 0) {
      // get, but need to compute since saved member is empty
      interfaces_ = InterfacesFromBoundaryCenters(SubPatchCenters(),
                                                  tolerance_,
                                                  ParaDim());
    }

    // regardless of set/get, it will always return the saved member of
    // current state
    return interfaces_;
  }

  /// @brief Gets IDs of boundary patch
  py::array_t<int> BoundaryPatchIds() {
    // return if it exists
    if (boundary_patch_ids_.size() != 0) {
      return boundary_patch_ids_;
    }

    // prepare only ingredient - interfaces
    // negative entries in interfaces mean boundary patch
    auto interfaces = Interfaces({});
    const int n_interfaces = interfaces.size();
    const int* interfaces_ptr = static_cast<int*>(interfaces.request().ptr);

    // create new array
    boundary_patch_ids_ = py::array_t<int>(n_interfaces);
    int* boundary_patch_ids_ptr =
        static_cast<int*>(boundary_patch_ids_.request().ptr);

    // add ids of negative interface entries
    int j{};
    for (int i{}; i < n_interfaces; ++i) {
      if (interfaces_ptr[i] < 0) {
        boundary_patch_ids_ptr[j] = i;
        ++j;
      }
    }

    // resize ids
    boundary_patch_ids_.resize({j}, false);

    return boundary_patch_ids_;
  }

  /// @brief Gets boundary multi patch
  std::shared_ptr<PyMultiPatch> BoundaryMultiPatch() {
    // return early if exists
    if (boundary_multi_patch_) {
      return boundary_multi_patch_;
    }

    // get boundary_patch_ids;
    auto boundary_pid = BoundaryPatchIds();
    int* boundary_pid_ptr = static_cast<int*>(boundary_pid.request().ptr);
    const int n_boundary_pid = boundary_pid.size();

    // values needed for offset
    const int n_subpatches = ParaDim() * 2;

    // create boundary patches with enough space
    CoreSplineVector boundary_core_patches(n_boundary_pid);

    // extract patches
    auto extract_boundaries_step = [&](const int begin, const int total_) {
      for (int i{begin}; i < total_; i += n_default_threads_) {
        const auto [i_spline, i_subpatch] =
            std::div(boundary_pid_ptr[i], n_subpatches);
        boundary_core_patches[i] =
            core_patches_[i_spline]->SplinepyExtractBoundary(i_subpatch);
      }
    };

    // Execute in parallel
    splinepy::utils::NThreadExecution(extract_boundaries_step,
                                      n_boundary_pid,
                                      n_default_threads_,
                                      splinepy::utils::NThreadQueryType::Step);

    // create return patch
    boundary_multi_patch_ = std::make_shared<PyMultiPatch>();
    boundary_multi_patch_->patches_ = ToPySplineList(boundary_core_patches);
    boundary_multi_patch_->core_patches_ = std::move(boundary_core_patches);

    return boundary_multi_patch_;
  }

  /// @brief returns derived class of boundary multi patch
  /// @return
  py::object PyBoundaryMultiPatch() { return ToDerived(BoundaryMultiPatch()); }

  /// @brief Evaluate at query points
  /// @param queries Query points
  /// @param nthreads Number of threads to use
  py::array_t<double> Evaluate(py::array_t<double> queries,
                               const int nthreads) {
    // use first spline as dimension guide line
    const int para_dim = ParaDim();
    const int dim = Dim();

    // query dim check
    CheckPyArrayShape(queries, {-1, para_dim}, true);

    // prepare input and output
    double* queries_ptr = static_cast<double*>(queries.request().ptr);
    const int n_splines = core_patches_.size();
    const int n_queries = queries.shape(0);
    const int n_total = n_splines * n_queries;
    py::array_t<double> evaluated({n_total, dim});
    double* evaluated_ptr = static_cast<double*>(evaluated.request().ptr);

    // each thread evaluates similar amount of queries from each spline
    auto evaluate_step = [&](int begin, int total_) {
      for (int i{begin}; i < total_; i += nthreads) {
        const auto [i_spline, i_query] = std::div(i, n_queries);
        core_patches_[i_spline]->SplinepyEvaluate(
            &queries_ptr[i_query * para_dim],
            &evaluated_ptr[(i_spline * n_queries + i_query) * dim]);
      }
    };

    // exe
    splinepy::utils::NThreadExecution(evaluate_step,
                                      n_total,
                                      nthreads,
                                      splinepy::utils::NThreadQueryType::Step);

    return evaluated;
  }

  /// @brief Sample multi patch
  /// @param resolution
  /// @param nthreads Number of threads to use
  /// @param same_parametric_bounds
  py::array_t<double> Sample(const int resolution,
                             const int nthreads,
                             const bool same_parametric_bounds) {
    const int para_dim = ParaDim();
    const int dim = Dim();

    // n_queries, and n_splines;
    int n_queries{1};
    const int n_splines = core_patches_.size();

    // prepare resolutions
    IntVector resolutions_vector(para_dim);
    int* resolutions = resolutions_vector.data();
    for (int i{}; i < para_dim; ++i) {
      resolutions[i] = resolution;
      n_queries *= resolution;
    }

    // prepare input /  output
    const int n_total = n_splines * n_queries;
    py::array_t<double> sampled({n_total, dim});
    double* sampled_ptr = static_cast<double*>(sampled.request().ptr);

    // if you know all the queries have same parametric bounds
    // you don't need to re-compute queries
    if (same_parametric_bounds) {

      // get para bounds
      DoubleVector para_bounds_vector(2 * para_dim);
      double* para_bounds = para_bounds_vector.data();
      core_patches_[0]->SplinepyParametricBounds(para_bounds);

      // prepare queries
      DoubleVector queries_vector(n_queries * para_dim);
      double* queries = queries_vector.data();

      // use grid point generator to fill queries
      splinepy::utils::CStyleArrayPointerGridPoints gp_generator(para_dim,
                                                                 para_bounds,
                                                                 resolutions);
      gp_generator.Fill(queries);

      // create lambda for nthread exe
      auto sample_same_bounds_step = [&](int begin, int total_) {
        for (int i{begin}; i < total_; i += nthreads) {
          const auto [i_spline, i_query] = std::div(i, n_queries);
          core_patches_[i_spline]->SplinepyEvaluate(
              &queries[i_query * para_dim],
              &sampled_ptr[(i_spline * n_queries + i_query) * dim]);
        }
      };

      splinepy::utils::NThreadExecution(
          sample_same_bounds_step,
          n_total,
          nthreads,
          splinepy::utils::NThreadQueryType::Step);

    } else {
      // here, we will execute 2 times:
      //   first, to create grid point helpers for each spline
      //   second, to sample

      // create a container to hold grid point helper.
      splinepy::utils::DefaultInitializationVector<
          splinepy::utils::CStyleArrayPointerGridPoints>
          grid_points(n_splines);

      // create grid_points
      auto create_grid_points = [&](int begin, int end) {
        DoubleVector para_bounds_vector(2 * para_dim);
        double* para_bounds = para_bounds_vector.data();

        for (int i{begin}; i < end; ++i) {
          // get para_bounds
          core_patches_[i]->SplinepyParametricBounds(para_bounds);
          // setup grid points helper
          grid_points[i].SetUp(para_dim, para_bounds, resolutions);
        }
      };

      // pre compute entries -> this one is a chunk query
      splinepy::utils::NThreadExecution(create_grid_points,
                                        n_splines,
                                        nthreads);

      // similar to the one with same_parametric_bounds, except it computes
      // query on the fly
      auto sample_step = [&](int begin, int total_) {
        // each thread needs just one query array
        DoubleVector thread_query_vector(para_dim);
        double* thread_query = thread_query_vector.data();

        for (int i{begin}; i < total_; i += nthreads) {
          const auto [i_spline, i_query] = std::div(i, n_queries);
          const auto& gp_helper = grid_points[i_spline];
          gp_helper.IdToGridPoint(i_query, thread_query);
          core_patches_[i_spline]->SplinepyEvaluate(
              thread_query,
              &sampled_ptr[(i_spline * n_queries + i_query) * dim]);
        }
      };

      // exe - this one is step
      splinepy::utils::NThreadExecution(
          sample_step,
          n_total,
          nthreads,
          splinepy::utils::NThreadQueryType::Step);
    }

    return sampled;
  }

  /// @brief Adds fields
  /// @param fields
  /// @param check_name
  /// @param check_dims
  /// @param check_degrees
  /// @param check_control_mesh_resolutions
  /// @param nthreads
  void AddFields(py::list& fields,
                 const bool check_name,
                 const bool check_dims,
                 const bool check_degrees,
                 const bool check_control_mesh_resolutions,
                 const int nthreads) {

    // allocate space
    const auto n_current_fields = static_cast<int>(field_multi_patches_.size());
    const auto n_new_fields = static_cast<int>(fields.size());
    const int n_total_fields = n_current_fields + n_new_fields;
    field_multi_patches_.resize(n_total_fields);

    // some hint values for null splines
    const int para_dim = ParaDim();
    const int dim = Dim();

    // prepare error message
    std::string field_mismatch_info{};
    std::mutex field_mismatch_mutex;

    // loop each field
    auto field_to_multi_patch = [&](int begin, int end) {
      // with in this nthread exe, we only use nthreads=1. This won't create
      // any threads.
      for (int i{n_current_fields + begin}, j{begin}; j < end; ++i, ++j) {
        // create multipatch
        py::list casted_list = fields[j].template cast<py::list>();
        auto field_multi_patch =
            std::make_shared<PyMultiPatch>(casted_list, para_dim, dim, 1);
        // propagate same_parametric_bounds_ flag
        field_multi_patch->same_parametric_bounds_ = same_parametric_bounds_;
        // set item
        field_multi_patches_[i] = field_multi_patch;

        try {
          // check mismatch - doesn't check null splines
          RaiseMismatch(core_patches_,
                        field_multi_patches_[i]->core_patches_,
                        check_name,
                        check_dims,
                        check_dims,
                        check_degrees,
                        check_control_mesh_resolutions,
                        1);

        } catch (const std::runtime_error& e) {
          // set true, and add error message.
          std::lock_guard<std::mutex> guard(field_mismatch_mutex);
          field_mismatch_info += "[mismatch error from the field with index ("
                                 + std::to_string(j) + ")]\n" + e.what();
        }
      }
    };

    // raise, if there were error.
    if (field_mismatch_info.size() != 0) {
      splinepy::utils::PrintAndThrowError(field_mismatch_info);
    }

    // all good, extend to list
    field_list_ += fields;
  }

  /// @brief Gets list of fields
  py::list GetFields() { return field_list_; }
};

/// @brief Add multi patch.
/// Returns [connectivity, vertex_ids, edge_information, boundaries]
/// @param m
inline void add_multi_patch(py::module& m) {
  m.def("interfaces_from_boundary_centers",
        &splinepy::py::InterfacesFromBoundaryCenters,
        py::arg("face_center_vertices"),
        py::arg("tolerance"),
        py::arg("para_dim"));
  m.def("extract_all_boundary_splines",
        &splinepy::py::ExtractAllBoundarySplines,
        py::arg("splines"),
        py::arg("interfaces"),
        py::arg("nthreads") = 1);
  m.def("orientations",
        &splinepy::py::GetBoundaryOrientations,
        py::arg("splines"),
        py::arg("base_ids"),
        py::arg("base_face_ids"),
        py::arg("neighbor_ids"),
        py::arg("neighbor_face_ids"),
        py::arg("tolerance"),
        py::arg("nthreads") = 1);
  m.def("boundaries_from_continuity",
        &splinepy::py::AddBoundariesFromContinuity,
        py::arg("boundary_splines"),
        py::arg("boundary_interfaces"),
        py::arg("global_interfaces"),
        py::arg("tolerance"),
        py::arg("nthreads") = 1);

  py::class_<splinepy::py::PyMultiPatch,
             std::shared_ptr<splinepy::py::PyMultiPatch>>
      klasse(m, "PyMultiPatch");

  klasse.def(py::init<>())
      .def(py::init<py::list&, const int, const bool>())
      .def(py::init<std::shared_ptr<PyMultiPatch>&>()) // for to_derived()
      .def_readwrite("n_default_threads", &PyMultiPatch::n_default_threads_)
      .def_readwrite("same_parametric_bounds",
                     &PyMultiPatch::same_parametric_bounds_)
      .def_readwrite("tolerance", &PyMultiPatch::tolerance_)
      .def("clear", &PyMultiPatch::Clear)
      .def_property_readonly("para_dim", &PyMultiPatch::ParaDim)
      .def_property_readonly("dim", &PyMultiPatch::Dim)
      .def_property_readonly("name", &PyMultiPatch::Name)
      .def_property_readonly("whatami", &PyMultiPatch::WhatAmI)
      .def_property("patches",
                    &PyMultiPatch::GetPatches,
                    &PyMultiPatch::SetPatchesDefault)
      .def("sub_multi_patch", &PyMultiPatch::PySubMultiPatch)
      .def("sub_patch_centers", &PyMultiPatch::SubPatchCenters)
      .def("interfaces", &PyMultiPatch::Interfaces)
      .def("boundary_patch_ids", &PyMultiPatch::BoundaryPatchIds)
      .def("boundary_multi_patch", &PyMultiPatch::PyBoundaryMultiPatch)
      .def("evaluate",
           &PyMultiPatch::Evaluate,
           py::arg("queries"),
           py::arg("nthreads"))
      .def("sample",
           &PyMultiPatch::Sample,
           py::arg("resolution"),
           py::arg("nthreads"),
           py::arg("same_parametric_bounds"))
      .def("add_fields",
           &PyMultiPatch::AddFields,
           py::arg("fields"),
           py::arg("check_name"),
           py::arg("check_dims"),
           py::arg("check_degrees"),
           py::arg("checK_control_mesh_resolutions"),
           py::arg("nthreads"))
      .def("fields", &PyMultiPatch::GetFields)
      //.def("", &PyMultiPatch::)
      ;
}

} // namespace splinepy::py
