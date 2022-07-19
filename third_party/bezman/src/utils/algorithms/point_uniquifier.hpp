/*
MIT License

Copyright (c) 2022 zwar@ilsb.tuwien.ac.at

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef UTILS_ALGORITHMS_POINT_UNIQUIFIER_HPP
#define UTILS_ALGORITHMS_POINT_UNIQUIFIER_HPP

#include <array>
#include <cassert>
#include <vector>

#include "bezman/src/bezier_spline_group.hpp"
#include "bezman/src/point.hpp"
#include "bezman/src/utils/algorithms/hypercube.hpp"
#include "bezman/src/utils/algorithms/sort.hpp"

namespace bezman::utils::algorithms {

/**
 * @brief Takes a set of points and determines point-connectivity using a metric
 *
 * Duplicate Points are not eliminated, assuming that aa maximum of two points
 * are equivalent. If this is not the case an exception is thrown. In theory
 * this has complexity O(nlogn) whereas a KDTree has complexity O(n (logn)^dim).
 *
 * @tparam physical_dimension   Entries in Point
 * @tparam ScalarType           Type of Entry in Point
 * @param face_center_points    List of intersection points (Mid-Point of
 *                              face-vertices)
 * @param orientation_metric    Vector along which the points are ordered
 *                              prior to their uniquification
 * @param tolerance             Tolerance for vector contraction
 * @return Connectivity         Element Connectivity
 */
template <std::size_t physical_dimension, typename ScalarType,
          std::size_t number_of_element_faces>
auto FindConnectivity(
    const std::vector<Point<physical_dimension, ScalarType>>&
        face_center_points,
    Point<physical_dimension, ScalarType> orientation_metric,
    const std::array<std::size_t, number_of_element_faces>& opposite_face_list,
    const ScalarType tolerance = 1e-5) {
  // Check if number of faces is a divisor of the point list length
  Logger::Logging("Determining connectivity by analyzing face centers");
  assert(("Wrong number of faces and center points",
          face_center_points.size() % number_of_element_faces == 0));

  // Assure Metric is normed and non-zero
  if (orientation_metric.EuclidianNorm() < 1e-20) {
    Logger::Warning(
        "Metric has no length. Chose non-zero "
        "metric for ordering points");
    Logger::Warning("Fall back to default metric, which is {1., 1., ...}");
    orientation_metric.fill(1.);
  }
  const Point<physical_dimension, ScalarType> normed_orientation_metric =
      orientation_metric *
      (static_cast<ScalarType>(1.) / orientation_metric.EuclidianNorm());

  // Store information in Auxiliary Values
  const std::size_t n_total_points{face_center_points.size()};
  const ScalarType tolerance_squared{tolerance * tolerance};

  // Init connectivity and metric value
  // (-1 : boundary, -2 : untouched)
  const std::size_t number_of_elements =
      face_center_points.size() / number_of_element_faces;
  std::vector<std::array<std::size_t, number_of_element_faces>> connectivity(
      // Size of vector
      number_of_elements,
      // Lambda function to initialize an array with constant value in size of
      // face-number
      []() {
        std::array<std::size_t, number_of_element_faces> a{};
        a.fill(static_cast<std::size_t>(-2));
        return a;
      }());
  std::vector<ScalarType> scalar_metric{};
  scalar_metric.reserve(n_total_points);

  // Check Metric Dimension and Vector Size
  for (unsigned int i{}; i < n_total_points; i++) {
    scalar_metric.push_back(normed_orientation_metric * face_center_points[i]);
  }

  // Sort Metric Vector
  const auto metric_order_indices = algorithms::IndexListSort(scalar_metric);

  // Loop over points
  for (unsigned int lower_limit = 0; lower_limit < n_total_points - 1;
       lower_limit++) {
    // Loop over all points regardless of whether they have been touched or not,
    // and then check the validity of the connection Point already processed
    bool found_duplicate = false;
    // Now check allowed range for duplicates
    unsigned int upper_limit = lower_limit + 1;
    while (scalar_metric[metric_order_indices[upper_limit]] -
                   scalar_metric[metric_order_indices[lower_limit]] <
               tolerance &&
           upper_limit < n_total_points) {
      // Check if the two points are duplicates
      found_duplicate = (face_center_points[metric_order_indices[lower_limit]] -
                         face_center_points[metric_order_indices[upper_limit]])
                            .SquaredEuclidianNorm() < tolerance_squared;
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
    // 2. The point it connects to must be on opposite sides on the neighboring
    //    element, else there is an orientation problem in the mesh and the mesh
    //    can not be exported in mfem format, this check needs to be disabled
    //    for general connectivity
    const std::size_t id_start_point = metric_order_indices[lower_limit];
    const std::size_t element_id_start =
        id_start_point / number_of_element_faces;
    const std::size_t element_face_id_start =
        id_start_point - element_id_start * number_of_element_faces;
    // is that id_start_point % number_of_element_faces?

    if (found_duplicate) {
      // Calculate indices
      const std::size_t id_end_point = metric_order_indices[upper_limit];
      const std::size_t element_id_end = id_end_point / number_of_element_faces;
      const std::size_t element_face_id_end =
          id_end_point - element_id_end * number_of_element_faces;
      // Check 1. (@todo EXCEPTION)
      assert(connectivity[element_id_start][element_face_id_start] ==
             static_cast<std::size_t>(-2));

      if (connectivity[element_id_start][element_face_id_start] !=
          static_cast<std::size_t>(-2)) {
        Logger::TerminatingError(
            "Connectivity connection is invalid. "
            "Found conflicting interceptions");
      }

      // Check 2. (@todo EXCEPTION)
      // TODO check if mfem format is used for the output -> if not do not check
      if (opposite_face_list[element_face_id_start] != element_face_id_end) {
        Logger::TerminatingError("Orientation Problem for MFEM-mesh output.");
      }
#ifdef NDEBUG
      if (opposite_face_list[element_face_id_end] != element_face_id_start) {
        Logger::Error("Orientation Problem for MFEM-mesh output.");
      }
#endif
      // If both tests passed, update connectivity
      connectivity[element_id_start][element_face_id_start] = element_id_end;
      connectivity[element_id_end][element_face_id_end] = element_id_start;
    } else {
      // set Boundary-ID
      if (connectivity[element_id_start][element_face_id_start] ==
          static_cast<std::size_t>(-2)) {
        connectivity[element_id_start][element_face_id_start] =
            static_cast<std::size_t>(-1);
      }
    }
  }

  // Treat last remaining point in scalar metric vector
  const std::size_t last_id = metric_order_indices[n_total_points - 1];
  const std::size_t last_element = last_id / number_of_element_faces;
  const std::size_t last_face =
      last_id - last_element * number_of_element_faces;
  if (connectivity[last_element][last_face] == static_cast<std::size_t>(-2)) {
    connectivity[last_element][last_face] = static_cast<std::size_t>(-1);
  }
  Logger::Logging("Found " + std::to_string(last_id) + " connections for " +
                  std::to_string(n_total_points) + " faces");
  return connectivity;
}

/**
 * @brief Finds duplicate Points and returns a list with indices that can be
 * used to build this list.
 *
 * example:
 * the list
 * [[0,1],[2,1],[0,1],[0,2]]
 * could return
 * [2,1,2,0] (order depends on orientation metric)
 */
template <std::size_t physical_dimension, typename ScalarType>
std::vector<std::size_t> IndexUniquePointList(
    const std::vector<Point<physical_dimension, ScalarType>>&
        original_point_list,
    const Point<physical_dimension, ScalarType> orientation_metric,
    const ScalarType tolerance = 1e-5) {
  Logger::Logging("Indexing unique point list");
  // Assure Metric is normed and non-zero
  if (orientation_metric.EuclidianNorm() <= 0) {
    Logger::TerminatingError("Metric is not normed or zero");
  }
  const Point<physical_dimension, ScalarType> normed_orientation_metric =
      orientation_metric *
      (static_cast<ScalarType>(1.) / orientation_metric.EuclidianNorm());

  // Store information in Auxiliary Values
  const std::size_t n_total_points{original_point_list.size()};
  const ScalarType tolerance_squared{tolerance * tolerance};

  // Init unique_indices and metric value
  // (-1 : untouched)
  // {in c++20 this expression could be constexpr}
  std::vector<std::size_t> unique_indices(
      // Size of vector
      n_total_points,
      // default value
      static_cast<std::size_t>(-1));

  // Initialize Metric
  std::vector<ScalarType> scalar_metric{};
  scalar_metric.reserve(n_total_points);

  // Check Metric Dimension and Vector Size
  for (unsigned int i{}; i < n_total_points; i++) {
    scalar_metric.push_back(normed_orientation_metric * original_point_list[i]);
  }

  // Sort Metric Vector
  const auto metric_order_indices = algorithms::IndexListSort(scalar_metric);

  // Start Uniquifying
  Logger::ExtendedInformation("Start unique indexing of control points");
  std::size_t number_of_new_points{};
  for (std::size_t lower_limit{0}; lower_limit < n_total_points - 1;
       lower_limit++) {
    // Point already processed
    if (unique_indices[metric_order_indices[lower_limit]] !=
        static_cast<std::size_t>(-1)) {
      continue;
    }

    // Point has not been processed -> add it to new point list
    unique_indices[metric_order_indices[lower_limit]] = number_of_new_points;

    // Now check allowed range for duplicates
    unsigned int upper_limit = lower_limit + 1;
    while (scalar_metric[metric_order_indices[upper_limit]] -
                   scalar_metric[metric_order_indices[lower_limit]] <
               tolerance &&
           upper_limit < n_total_points) {
      const bool found_duplicate =
          (original_point_list[metric_order_indices[lower_limit]] -
           original_point_list[metric_order_indices[upper_limit]])
              .SquaredEuclidianNorm() < tolerance_squared;
      if (found_duplicate) {
        if (unique_indices[metric_order_indices[upper_limit]] !=
            static_cast<std::size_t>(-1)) {
          Logger::TerminatingError(
              "Failure in indexing Unique Point List. "
              "Found more than two different indices in less than two "
              "tolerances proximity");
        }
        unique_indices[metric_order_indices[upper_limit]] =
            number_of_new_points;
      }
      upper_limit++;
    }
    number_of_new_points++;
  }

  // Special case
  const auto& last_index = metric_order_indices.size() - 1;
  if (unique_indices[metric_order_indices[last_index]] ==
      static_cast<std::size_t>(-1)) {
    unique_indices[metric_order_indices[last_index]] = number_of_new_points;
  }
  Logger::Logging("Found " + std::to_string(number_of_new_points) +
                  " unique points out of " + std::to_string(n_total_points) +
                  " points");
  return unique_indices;
}

/**
 * @brief Get the Connectivity between splines in a SplineGroup
 *
 * Calculates the face points and uses FindConnectivity to identify
 * neighbors
 */
template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
auto GetConnectivityForSplineGroup(
    const BezierSplineGroup<parametric_dimension, PhysicalPointType,
                            ScalarType>& spline_group) {
  Logger::Logging("Determining connectivity");
  // Current implementation is only made for bi- and trivariates
  static_assert((parametric_dimension == 3 || parametric_dimension == 2),
                "High-Dimensional and Line Patches are not supported");

  // Array that stores opposite faces
  constexpr auto opposite_faces =
      HyperCube<parametric_dimension>::GetOppositeFaces();

  // Create Face-Center-Point Vector
  std::vector<PhysicalPointType> face_edges(spline_group.size() *
                                            opposite_faces.size());

  // Start Loop
  // (Instead of using the mean of the face vertices, using the sum)
  const std::size_t number_of_splines = spline_group.size();
  const std::size_t number_of_element_faces = opposite_faces.size();

  // Retrieve SubElementFace-Vertex Ids in local system to start calculating
  // face-mid-point
  constexpr auto subelement_vertex_ids =
      HyperCube<parametric_dimension>::SubElementVerticesToFace();

  for (std::size_t i_spline{}; i_spline < number_of_splines; i_spline++) {
    const auto global_vertex_id =
        HyperCube<parametric_dimension>::VertexIdForDegrees(
            spline_group[i_spline].GetDegrees());
    for (std::size_t i_face{}; i_face < number_of_element_faces; i_face++) {
      face_edges[i_spline * number_of_element_faces + i_face] =
          spline_group[i_spline].control_points
              [global_vertex_id[subelement_vertex_ids[i_face][0]]];
      for (std::size_t i_point{1}; i_point < subelement_vertex_ids[0].size();
           i_point++) {
        face_edges[i_spline * number_of_element_faces + i_face] +=
            spline_group[i_spline].control_points
                [global_vertex_id[subelement_vertex_ids[i_face][i_point]]];
      }
    }
  }

  // Get Connectivity
  return FindConnectivity(
      face_edges,
      // Metric for internal ordering
      spline_group.MaximumCorner() - spline_group.MinimumCorner(),
      opposite_faces);
}

}  // namespace bezman::utils::algorithms
#endif  // UTILS_UNIQUIFY_POINT_UNIQUIFIER_HPP
