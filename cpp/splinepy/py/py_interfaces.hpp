#pragma once

// pybind11
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "splinepy/utils/print.hpp"

namespace splinepy::py {

namespace py = pybind11;

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
  const double tolerance_squared{tolerance * tolerance};
  const int number_of_patches =
      face_center_vertices.size() / number_of_element_faces;
  const int physical_dimension = face_center_vertices.shape()[1];
  const int number_of_center_vertices{
      static_cast<int>(face_center_vertices.size())};
  const double* metric_ptr = static_cast<double*>(metric.request().ptr);
  const double* face_center_vertices_ptr =
      static_cast<double*>(face_center_vertices.request().ptr);

  // Consistency check
  if (!(face_center_vertices.shape()[0] % number_of_element_faces == 0)) {
    splinepy::utils::PrintAndThrowError(
        "Number of corner vertices invalid. Must be a multiple of the number "
        "of vertices per patch");
  }
  if (!(metric.size() == face_center_vertices.shape()[1])) {
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
  py::array_t<int> connectivity(number_of_patches * parametric_dimension * 2);
  int* connectivity_ptr = static_cast<int*>(connectivity.request().ptr);
  // Init connectivity and metric value
  // (-1 : boundary, -2 : untouched)
  std::fill(connectivity_ptr,
            connectivity_ptr + number_of_patches * parametric_dimension * 2,
            -2);

  // Check for Wrong number of faces and center points
  assert(number_of_center_vertices % number_of_element_faces == 0);

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
    double squared_euclidian_distance{};
    for (int i_phys{}; i_phys < physical_dimension; i_phys++) {
      const double distance_c =
          face_center_vertices_ptr[i_end * physical_dimension + i_phys]
          - face_center_vertices_ptr[i_start * physical_dimension + i_phys];
      squared_euclidian_distance += distance_c * distance_c;
    }
    return squared_euclidian_distance;
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
           && scalar_metric[metric_order_indices[upper_limit]]
                      - scalar_metric[metric_order_indices[lower_limit]]
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
      if (connectivity_ptr[element_id_start * number_of_element_faces
                           + element_face_id_start]
          != -2) {
        splinepy::utils::PrintAndThrowError(
            "Found conflicting interceptions, where more than two points are "
            "in the same position");
      }

      // If both tests passed, update connectivity
      connectivity_ptr[element_id_start * number_of_element_faces
                       + element_face_id_start] = element_id_end;
      connectivity_ptr[element_id_end * number_of_element_faces
                       + element_face_id_end] = element_id_start;
    } else {
      // set Boundary-ID
      if (connectivity_ptr[element_id_start * number_of_element_faces
                           + element_face_id_start]
          == -2) {
        connectivity_ptr[element_id_start * number_of_element_faces
                         + element_face_id_start] = -1;
      }
    }
  }

  // Treat last remaining point in scalar metric vector
  const std::size_t last_id =
      metric_order_indices[number_of_center_vertices - 1];
  const std::size_t last_element = last_id / number_of_element_faces;
  const std::size_t last_face =
      last_id - last_element * number_of_element_faces;
  if (connectivity_ptr[last_element * number_of_element_faces + last_face]
      == -2) {
    connectivity_ptr[last_element * number_of_element_faces + last_face] = -1;
  }

  // Resize buffer
  connectivity.resize({number_of_patches, number_of_element_faces});
  return connectivity;
}

// Provide function to add to module
inline void add_interfaces(py::module& m) {
  m.def("interfaces_from_boundary_centers",
        &splinepy::py::FindConnectivityFromCenters,
        py::arg("boundary_centers"),
        py::arg("parametric_dimension"),
        py::arg("metric"),
        py::arg("tolerance"));
}

} // namespace splinepy::py
