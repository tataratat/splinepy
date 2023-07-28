#pragma once

#include <array>
#include <type_traits>

#include "splinepy/splines/bezier.hpp"
#include "splinepy/splines/rational_bezier.hpp"
#include "splinepy/splines/splinepy_base.hpp"

namespace splinepy::splines::helpers {

/**
 * @brief Extract Bezier Patch IDs
 *
 * This function assumed that repeated knots have been inserted, such that the
 * continuity at each knot is C0. It returns the indices that are associated to
 * the individual bezier patches and returns a vector of ID lists that can be
 * used to create the new splines
 *
 * @tparam IndexingType integral type that represents the IDs - does not change
 * return type
 * @param degrees pointer to degrees vector
 * @param n_patches_per_para_dim pointer to vector with same dimension that
 * represents the number of patches within the B-Spline type along a respective
 * parametric dimension
 * @param para_dim parametric dimension
 * @return std::vector<std::vector<int>>
 */
template<typename IndexingType>
std::vector<std::vector<int>>
ExtractBezierPatchIDs(const IndexingType* degrees,
                      const int* n_patches_per_para_dim,
                      const int& para_dim) {
  // Check for type
  static_assert(std::is_integral_v<IndexingType>,
                "Unsupported type for ExtractBezierPatchIDs");
  // Number of total patches and ctps per patch
  int n_total_patches = n_patches_per_para_dim[0];
  IndexingType n_ctps_per_patch = degrees[0] + 1;

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
    for (int i_local_id{}; i_local_id < static_cast<int>(n_ctps_per_patch);
         i_local_id++) {
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

/**
 * @brief Extracts Bezier patches of a B-Spline/NURBS type
 *
 * @tparam as_base flag to determine input type avoid circular dependency
 * @tparam SplineType Spline-type (NURBS of BSpline)
 * @param input spline to be separated
 * @return auto vector of Bezier types, either Rational or polynomial
 */
template<bool as_base = false, typename SplineType>
auto ExtractBezierPatches(SplineType& input) {
  // Start by identifying types
  constexpr int para_dim = SplineType::kParaDim;
  constexpr int dim = SplineType::kDim;
  constexpr bool is_rational = SplineType::kIsRational;

  using ReturnType =
      std::conditional_t<is_rational,
                         splinepy::splines::RationalBezier<para_dim, dim>,
                         splinepy::splines::Bezier<para_dim, dim>>;
  using ReturnVectorValueType =
      std::conditional_t<as_base, splinepy::splines::SplinepyBase, ReturnType>;
  using PointType = typename ReturnType::Coordinate_;

  // Predetermine some auxiliary values
  // Access spline degrees and determine offsets
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
            bsplinelib::Dimension{i_p_dim});
    n_patches_per_para_dim[i_p_dim] = std::get<0>(bezier_information);
    n_ctps_per_para_dim[i_p_dim] =
        n_patches_per_para_dim[i_p_dim] * degrees[i_p_dim] + 1;
    const auto& knot_vector_ref = std::get<1>(bezier_information);

    // Insert knot into the copy of the spline before extraction
    // this is the most costly part of the calculation
    for (std::size_t i_knot{}; i_knot < knot_vector_ref.size(); i_knot++) {
      input.InsertKnot(bsplinelib::Dimension{i_p_dim}, knot_vector_ref[i_knot]);
    }
  }
  // Auxiliary function
  const auto& ControlPointVector = [&](const int& id) {
    if constexpr (is_rational) {
      return input.GetWeightedVectorSpace()[bsplinelib::Index{id}];
    } else {
      return input.GetVectorSpace()[bsplinelib::Index{id}];
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
  std::vector<std::shared_ptr<ReturnVectorValueType>> bezier_list{};
  bezier_list.reserve(n_total_patches);

  // Loop over the individual patches
  for (int i_patch{}; i_patch < n_total_patches; i_patch++) {

    // Init vectors required for initialization
    std::vector<PointType> ctps;
    ctps.reserve(n_ctps_per_patch);
    std::vector<double> weights;
    if constexpr (is_rational) {
      weights.reserve(n_ctps_per_patch);
    }
    const auto& local_ids = ids[i_patch];

    // Extract relevant coordinates
    for (int i_local_id{}; i_local_id < n_ctps_per_patch; i_local_id++) {
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
      if constexpr (is_rational) {
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
    if constexpr (is_rational) {
      bezier_list.push_back(std::make_shared<ReturnType>(
          // degrees
          degrees,
          // CTPS non-weighted
          ctps,
          // weights
          weights));
    } else {
      bezier_list.push_back(std::make_shared<ReturnType>(
          // degrees
          degrees,
          // ctps
          ctps));
    }
  }
  return bezier_list;
}

} /* namespace splinepy::splines::helpers */
