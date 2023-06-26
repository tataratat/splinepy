#pragma once

#include <array>
#include <type_traits>

#include "bezman/src/bezier_spline.hpp"
#include "bezman/src/point.hpp"
#include "bezman/src/rational_bezier_spline.hpp"
#include "splinepy/splines/bspline.hpp"
#include "splinepy/splines/type_identifier.hpp"

namespace splinepy::splines {

/// Fills raw pointer with degrees. applicable for both nurbs and bspline
template<typename SplineType, typename DegreeType>
inline void FillDegrees(const SplineType& spline, DegreeType* degree_ptr) {
  const auto& parameter_space = spline.GetParameterSpace();
  for (int i{}; i < SplineType::kParaDim; ++i) {
    degree_ptr[i] = parameter_space.GetDegrees()[i].Get();
  }
}

/// Fills raw pointer with knot_vector. applicable for both nurbs and bspline
template<typename SplineType, typename IndexT, typename KnotType>
inline void FillKnotVector(const SplineType& spline,
                           IndexT para_dim,
                           KnotType* knot_vector) {
  const auto& parameter_space = spline.GetParameterSpace();
  const auto& requested_knot_vector =
      *parameter_space.GetKnotVectors()[para_dim];

  for (int i{}; i < requested_knot_vector.GetSize(); ++i) {
    knot_vector[i] = requested_knot_vector[splinelib::Index{i}].Get();
  }
}

/// Returns number of control_points
template<typename SplineType>
inline int GetNumberOfControlPoints(const SplineType& spline) {
  return spline.GetVectorSpace().GetNumberOfCoordinates();
}

// Naturally here would be a place for FillControlPoints.
// However, bspline and nurbs need different approaches for that.
// Hence, not here.

/// Returns parametric bounds as (para_dim x 2) array.
template<typename SplineType>
inline std::array<std::array<double, SplineType::kParaDim>, 2>
GetParametricBounds(const SplineType& spline) {
  std::array<std::array<double, SplineType::kParaDim>, 2> parametric_bounds{};

  const auto& parameter_space = spline.GetParameterSpace();
  size_t i{};
  for (const auto& knotvector : parameter_space.GetKnotVectors()) {
    parametric_bounds[0][i] = knotvector->GetFront().Get();
    parametric_bounds[1][i] = knotvector->GetBack().Get();
    ++i;
  }

  return parametric_bounds;
}

std::vector<std::vector<int>>
ExtractBezierPatchIDs(const int* degrees,
                      const int* n_patches_per_para_dim,
                      const int& para_dim) {
  // Number of total patches and ctps per patch
  int n_total_patches = n_patches_per_para_dim[0];
  int n_ctps_per_patch = degrees[0] + 1;

  // Offsets for start values of individual patches
  std::vector<std::size_t> bezier_index_offsets{};
  bezier_index_offsets.reserve(para_dim);
  bezier_index_offsets.push_back(1);
  for (int i_para_dim{1}; i_para_dim < para_dim; i_para_dim++) {
    n_total_patches *= n_patches_per_para_dim[i_para_dim];
    n_ctps_per_patch *= degrees[i_para_dim] + 1;
    bezier_index_offsets.push_back(bezier_index_offsets[i_para_dim - 1]
                                   * (degrees[i_para_dim - 1] + 1));
  }

  // Init return values
  std::vector<std::vector<int>> list_of_id_lists(n_total_patches);
  for (int i_patch{}; i_patch < n_total_patches; i_patch++) {
    // Determine internal postitions in local coord system
    std::vector<std::size_t> patch_ctp_id_offsets{};
    patch_ctp_id_offsets.reserve(para_dim);
    int ii{i_patch};
    // Determine the parameter wise ids of the patch (i.e. the
    // patch-coordinates) and calculate the required index offsets
    for (int i{}; i < para_dim; i++) {
      // ID in spline coordinate system of current patch
      const int patch_coord = static_cast<int>(ii % n_patches_per_para_dim[i]);
      ii -= patch_coord;
      ii /= n_patches_per_para_dim[i];
      // Coordinate offset of the control points indices
      patch_ctp_id_offsets.push_back(patch_coord * degrees[i]);
    }

    // Init vectors required for initialization
    std::vector<int>& ids = list_of_id_lists[i_patch];
    ids.reserve(n_ctps_per_patch);

    // Extract relevant coordinates
    for (int i_local_id{}; i_local_id < n_ctps_per_patch; i_local_id++) {
      int global_id{};
      int n_ctps_in_previous_layers{1};
      // Determine index of local point in global spline
      for (int i_para_dim{}; i_para_dim < para_dim; i_para_dim++) {
        // First id in local system
        const int local_id = (i_local_id / bezier_index_offsets[i_para_dim])
                             % (degrees[i_para_dim] + 1);
        // Add patch offsets
        global_id += (local_id + patch_ctp_id_offsets[i_para_dim])
                     * n_ctps_in_previous_layers;
        // Multiply to index offset
        n_ctps_in_previous_layers *=
            n_patches_per_para_dim[i_para_dim] * degrees[i_para_dim] + 1;
      }
      ids.push_back(global_id);
    }
  }
  return list_of_id_lists;
}

template<typename SplineType,
         std::enable_if_t<isBSpline_v<SplineType> || isNurbs_v<SplineType>>* =
             nullptr>
auto ExtractBezierPatches(SplineType& input) {
  // Start by identifying types
  constexpr int para_dim = SplineType::para_dim_;
  constexpr int dim = SplineType::dim_;
  using PointType =
      std::conditional_t<dim == 1, double, bezman::Point<dim, double>>;
  using ReturnType = std::conditional_t<
      isBSpline_v<SplineType>,
      bezman::BezierSpline<para_dim, PointType, double>,
      bezman::RationalBezierSpline<para_dim, PointType, double>>;

  // Predetermine some auxiliary values
  // Access Spline Degres and determine offsets
  std::array<std::size_t, para_dim> degrees{};
  for (int i_p{}; i_p < para_dim; i_p++) {
    degrees[i_p] = static_cast<std::size_t>(
        input.GetParameterSpace().GetDegrees()[i_p].Get());
  }
  std::array<int, para_dim> n_patches_per_para_dim{};
  std::array<int, para_dim> n_ctps_per_para_dim{};

  // Identify all internal knots and the number of required Bezier patches
  const auto parameter_space = input.GetParameterSpace();
  for (int i_p_dim{}; i_p_dim < para_dim; i_p_dim++) {
    // Extract internal knots
    const auto bezier_information =
        parameter_space.DetermineBezierExtractionKnots(
            // Use SplineLib Type
            splinelib::Dimension{i_p_dim});
    n_patches_per_para_dim[i_p_dim] = std::get<0>(bezier_information);
    n_ctps_per_para_dim[i_p_dim] =
        n_patches_per_para_dim[i_p_dim] * degrees[i_p_dim] + 1;
    const auto& knot_vector_ref = std::get<1>(bezier_information);

    // Insert knot into the copy of the spline before extraction
    // this is the most costly part of the calculation
    for (std::size_t i_knot{}; i_knot < knot_vector_ref.size(); i_knot++) {
      input.InsertKnot(splinelib::Dimension{i_p_dim}, knot_vector_ref[i_knot]);
    }
  }
  // Auxiliary function
  const auto& ControlPointVector = [&](const int& id) {
    if constexpr (isNurbs_v<SplineType>) {
      return input.GetWeightedVectorSpace()[splinelib::Index{id}];
    } else {
      return input.GetVectorSpace()[splinelib::Index{id}];
    }
  };

  // Retrieve id information
  const std::vector<std::vector<int>>& ids =
      ExtractBezierPatchIDs(degrees.data(),
                            n_patches_per_para_dim.data(),
                            para_dim);

  // Number of total patches
  const int n_total_patches = ids.size();
  const int n_ctps_per_patch = ids[0].size();

  // Init return value
  std::vector<ReturnType> bezier_list{};
  bezier_list.reserve(n_total_patches);

  // Loop over the individual patches
  for (int i_patch{}; i_patch < n_total_patches; i_patch++) {

    // Init vectors required for initialization
    std::vector<PointType> ctps;
    ctps.reserve(n_ctps_per_patch);
    std::vector<double> weights;
    if constexpr (isNurbs_v<SplineType>) {
      weights.reserve(n_ctps_per_patch);
    }
    const auto& local_ids = ids[i_patch];

    // Extract relevant coordinates
    for (std::size_t i_local_id{}; i_local_id < n_ctps_per_patch;
         i_local_id++) {
      const std::size_t global_id = local_ids[i_local_id];
      // Retrieve Control Point
      PointType point{};
      // Check if scalar (double instead of array)
      if constexpr (dim == 1) {
        point = ControlPointVector(global_id)[0];
      } else {
        for (int i_dim{}; i_dim < dim; i_dim++) {
          point[i_dim] = ControlPointVector(global_id)[i_dim];
        }
      }
      if constexpr (isNurbs_v<SplineType>) {
        const double weight = ControlPointVector(global_id)[dim];
        weights.push_back(weight);
        // control points are weighted in non-rational splines to facilitate
        // calculations. They therefore need to be divided by the weight to get
        // the true control point positions
        point /= weight;
      }
      ctps.push_back(point);
    }

    // Create Spline and add it to list
    if constexpr (isNurbs_v<SplineType>) {
      bezier_list.push_back(ReturnType(
          // degrees
          degrees,
          // CTPS non-weighted
          ctps,
          // weights
          weights));
    } else {
      bezier_list.push_back(ReturnType(
          // degrees
          degrees,
          // ctps
          ctps));
    }
  }
  return bezier_list;
}
} // namespace splinepy::splines
