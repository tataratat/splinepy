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
std::vector<std::vector<int>> ExtractBezierPatchIDs(
    const std::shared_ptr<splinepy::splines::SplinepyBase>& spline,
    const IndexingType* degrees,
    const int* n_patches_per_para_dim) {
  // Check for type
  static_assert(std::is_integral_v<IndexingType>,
                "Unsupported type for ExtractBezierPatchIDs");

  // use runtime call to support this function in py_knot_insertion_matrix
  if (!spline->SplinepyHasKnotVectors()) {
    splinepy::utils::PrintAndThrowError(
        spline->SplinepyWhatAmI(),
        "is not a valid type for ExtractBezierPatchIDs");
  }
  const int para_dim = spline->SplinepyParaDim();

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
  std::vector<std::size_t> patch_ctp_id_offsets(para_dim, 0);
  std::vector<std::vector<int>> list_of_id_lists(n_total_patches);
  const std::vector<std::vector<int>> multiplicities =
      spline->SplinepyKnotMultiplicities();
  for (int i_patch{}; i_patch < n_total_patches; i_patch++) {
    // Determine internal positions in local coord system
    int ii{i_patch};
    // Determine the parameter wise ids of the patch (i.e. the
    // patch-coordinates) and calculate the required index offsets
    for (int i{}; i < para_dim; i++) {
      // ID in spline coordinate system of current patch
      const int patch_coord = static_cast<int>(ii % n_patches_per_para_dim[i]);

      // Coordinate offset of the control points indices
      patch_ctp_id_offsets[i] =
          (patch_coord == 0)
              ? 0
              : patch_ctp_id_offsets[i] + multiplicities[i][patch_coord];

      ii -= patch_coord;
      ii /= n_patches_per_para_dim[i];
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
template<typename SplineType>
std::vector<std::shared_ptr<splinepy::splines::SplinepyBase>>
ExtractBezierPatches(const SplineType& spline) {
  static_assert(SplineType::kHasKnotVectors,
                "Invalid type for ExtractBezierPatches");
  // make copy of input spline
  auto input_ptr = spline.SplinepyDeepCopy();
  SplineType& input = *(std::dynamic_pointer_cast<SplineType>(input_ptr));

  // Start by identifying types
  constexpr int para_dim = SplineType::kParaDim;
  constexpr bool is_rational = SplineType::kIsRational;
  const int dim = input.SplinepyDim();

  // Predetermine some auxiliary values
  // Access spline degrees and determine offsets
  std::array<int, para_dim> degrees{};
  input.SplinepyCurrentProperties(degrees.data(), nullptr, nullptr, nullptr);

  std::array<int, para_dim> n_patches_per_para_dim{};
  std::array<int, para_dim> n_ctps_per_para_dim{};

  // Identify all internal knots and the number of required Bezier patches
  const auto& parameter_space = input.GetParameterSpace();
  for (int i_p_dim{}; i_p_dim < para_dim; i_p_dim++) {
    // Need to use BSplineLib's infamous NamedType
    const bsplinelib::Dimension pd_query{i_p_dim};

    // Extract internal knots
    const auto bezier_information =
        parameter_space.DetermineBezierExtractionKnots(pd_query);
    n_patches_per_para_dim[i_p_dim] = std::get<0>(bezier_information);
    n_ctps_per_para_dim[i_p_dim] =
        n_patches_per_para_dim[i_p_dim] * degrees[i_p_dim] + 1;
    const auto& required_knot_insertions = std::get<1>(bezier_information);

    // Insert knot into the copy of the spline before extraction
    // this is the most costly part of the calculation
    for (const auto& knot_to_insert : required_knot_insertions) {
      input.InsertKnot(pd_query, knot_to_insert);
    }
  }

  // Retrieve id information
  const std::vector<std::vector<int>>& ids =
      ExtractBezierPatchIDs(input_ptr,
                            degrees.data(),
                            n_patches_per_para_dim.data());

  // Number of total patches
  const int n_total_patches = ids.size();
  const int n_ctps_per_patch = ids[0].size();

  // Init return value
  std::vector<std::shared_ptr<splinepy::splines::SplinepyBase>> bezier_list{};
  bezier_list.reserve(n_total_patches);

  // get base control points to copy - in case of NURBS, this will be weighted
  // coordinates with the weights at the last column
  const auto& bspline_cps = input.GetCoordinates();

  // prepare temporary space to copy control point values
  splinepy::utils::Array<double, 2> extracted_cps(n_ctps_per_patch, dim);
  splinepy::utils::Array<double> extracted_weights;
  double* extracted_weights_data{nullptr};
  if constexpr (is_rational) {
    extracted_weights.Reallocate(n_ctps_per_patch);
    extracted_weights_data = extracted_weights.data();
  }

  // Loop over the individual patches
  for (const std::vector<int>& cp_ids : ids) {

    // loop over cp ids to extract
    int local_counter{};
    for (const int& cp_id : cp_ids) {
      if constexpr (is_rational) {
        // we have to un-weight
        // get beginning of from and to row
        const double* from_ptr = &bspline_cps(cp_id, 0);
        double* to_ptr = &extracted_cps(local_counter, 0);
        // invert weight
        const double w_inv = 1. / from_ptr[dim];

        // apply inverted weight and save to temporary
        for (int i{}; i < dim; ++i) {
          to_ptr[i] = from_ptr[i] * w_inv;
        }
        // copy weight
        extracted_weights[local_counter] = from_ptr[dim];

        ++local_counter;
      } else {
        // for non-rational bspline, we can just copy values
        std::copy_n(&bspline_cps(cp_id, 0),
                    dim,
                    &extracted_cps(local_counter++, 0));
      }
    }

    // use dynamic spline creator to create beziers.
    // weights's data will be nullptr if it is non rational.
    // all the values are copied by bezier and rational bezier
    bezier_list.push_back(splinepy::splines::SplinepyBase::SplinepyCreate(
        para_dim,
        dim,
        degrees.data(),
        nullptr,
        extracted_cps.data(),
        extracted_weights_data));
  }

  return bezier_list;
}

} /* namespace splinepy::splines::helpers */
