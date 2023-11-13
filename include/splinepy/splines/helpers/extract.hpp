#pragma once

#include <array>
#include <vector>

#include "splinepy/splines/helpers/properties.hpp"
#include "splinepy/splines/splinepy_base.hpp"
#include "splinepy/utils/arrays.hpp"
#include "splinepy/utils/grid_points.hpp"
#include "splinepy/utils/print.hpp"

namespace splinepy::splines::helpers {

/// Determine boundary ID from axis and extrem-value
inline int ExtractBoundaryFromAxisAndExtrema(const int& axis,
                                             const int& extreme) {
  // Determine corresponding ID
  return (extreme > 0) ? 2 * axis + 1 : 2 * axis;
}

template<typename SplineType, typename VectorType>
std::shared_ptr<splinepy::splines::SplinepyBase>
ExtractControlMeshSliceFromIDs(const SplineType& spline,
                               const VectorType& indices,
                               const int& plane_normal_axis) {
  // Perform sanity check
  if constexpr (SplineType::kParaDim == 1) {
    splinepy::utils::PrintWarning(
        "Sorry, we don't support control mesh slicing"
        "of 1-Parametric Dim splines. Returning empty spline.");
    return std::make_shared<SplineType>();
  } else {
    // boundary spline type
    using SelfBoundary = typename SplineType::BoundaryType_;

    // prepare return
    std::shared_ptr<SelfBoundary> boundary_spline;

    // create new one
    if constexpr (SplineType::kHasKnotVectors) {
      // hola BSpline families
      using VSpace = typename SelfBoundary::PhysicalSpace_;

      // slices splines can use parent's parameter space without
      // recreating basis.
      // RemoveOneParametricDimension() creates a new parameter space with
      // same shared pointers to all the required properties.
      auto pspace = spline.GetParameterSpace().RemoveOneParametricDimension(
          plane_normal_axis);

      // get size - dim should come from the coordinate itself, since nurbs have
      // dim + 1 shape
      const auto& from_coordinates = spline.GetCoordinates();
      const int dim = from_coordinates.Shape()[1];
      const int n_cps = indices.size();

      // VSpace - maybe it is also
      // create (weighted) vector space and allocate space for cps
      auto vspace = std::make_shared<VSpace>();
      auto& coords = vspace->GetCoordinates();
      coords.Reallocate(n_cps * dim);
      coords.SetShape(n_cps, dim);

      // copy loop - copy each rows
      int local_counter{};
      for (const auto& id : indices) {
        std::copy_n(&from_coordinates(id, 0), dim, &coords(local_counter++, 0));
      }

      // assign boundary spline - uses base' ctor
      boundary_spline = std::make_shared<SelfBoundary>(pspace, vspace);
    } else {
      // hola bezier families
      // form degrees
      using Degrees = typename SelfBoundary::Degrees_;
      Degrees b_degrees{};
      std::size_t ncps{1};
      int counter{};
      for (int i{}; i < SplineType::kParaDim; ++i) {
        if (i == plane_normal_axis) {
          continue;
        }
        const auto& ds = spline.GetDegrees();
        const auto& d = ds[i];
        b_degrees[counter] = d;
        ncps *= static_cast<std::size_t>(d + 1);
        counter++;
      }

      if constexpr (SplineType::kIsRational) {
        // only way to avoid unweight -> reweight is to init first with degrees
        // assign boundary spline
        boundary_spline = std::make_shared<SelfBoundary>(b_degrees);
        auto& b_cps = boundary_spline->GetWeightedControlPoints();
        auto& b_ws = boundary_spline->GetWeights();

        // copy
        for (std::size_t i{}; i < ncps; ++i) {
          const auto& id = indices[i];
          b_cps[i] = spline.GetWeightedControlPoints()[id];
          b_ws[i] = spline.GetWeights()[id];
        }
      } else {
        // non-rational bezier
        using Coordinates = typename SelfBoundary::Coordinates_;
        Coordinates b_coords;
        b_coords.reserve(ncps);
        for (const auto& id : indices) {
          // copy control points
          b_coords.push_back(spline.GetControlPoints()[id]);
        }

        // assign boundary spline
        boundary_spline = std::make_shared<SelfBoundary>(b_degrees, b_coords);
      }
    }
    // return here, this way compiler is happy
    return boundary_spline;
  }
}

template<typename SplineType>
std::shared_ptr<splinepy::splines::SplinepyBase>
ExtractBoundaryMeshSlice(const SplineType& spline, const int& boundary_id) {
  if constexpr (SplineType::kParaDim == 1) {
    splinepy::utils::PrintWarning(
        "Sorry, we don't support control mesh slicing"
        "of 1-Parametric Dim splines. Returning empty spline.");
    return std::make_shared<SplineType>();
  } else {

    // get ids on boundary
    const auto cmr =
        splinepy::splines::helpers::template GetControlMeshResolutions<int>(
            spline);
    const int plane_normal_axis = static_cast<int>(boundary_id / 2);
    const int plane_id =
        ((boundary_id % 2) == 0) ? 0 : cmr[plane_normal_axis] - 1;
    const auto ids_on_boundary = splinepy::utils::GridPoints::IdsOnHyperPlane(
        cmr.data(),
        static_cast<int>(cmr.size()),
        plane_normal_axis,
        plane_id);
    return ExtractControlMeshSliceFromIDs(spline,
                                          ids_on_boundary,
                                          plane_normal_axis);
  }
}

/// returns boundary spline, which has one less para_dim.
template<typename SplineType>
std::shared_ptr<splinepy::splines::SplinepyBase>
ExtractControlMeshSlice(const SplineType& spline,
                        const int& plane_normal_axis,
                        const int& plane_id) {
  if constexpr (SplineType::kParaDim == 1) {
    splinepy::utils::PrintWarning(
        "Sorry, we don't support control mesh slicing"
        "of 1-Parametric Dim splines. Returning empty spline.");
    return std::shared_ptr<
        typename SplineType::template SelfTemplate_<SplineType::kParaDim,
                                                    SplineType::kDim>>{};
  } else {
    // get ids on boundary
    const auto cmr =
        splinepy::splines::helpers::GetControlMeshResolutions(spline);
    const auto ids_on_boundary = splinepy::utils::GridPoints::IdsOnHyperPlane(
        cmr.data(),
        static_cast<int>(cmr.size()),
        plane_normal_axis,
        plane_id);
    return ExtractControlMeshSliceFromIDs(spline,
                                          ids_on_boundary,
                                          plane_normal_axis);
  }
}

/**
 * @brief Extract Bezier Patch IDs
 *
 * This function assumed that repeated knots have been inserted, such that the
 * continuity at each knot is C0. It returns the indices that are associated to
 * the individual bezier patches and returns a vector of ID lists that can be
 * used to create the new splines
 *
 * @param knot_multiplicities multiplicies of each unique knot per dim
 * @param degrees pointer to degrees vector
 * @param n_patches_per_para_dim pointer to vector with same dimension that
 * represents the number of patches within the B-Spline type along a respective
 * parametric dimension
 * @return std::vector<std::vector<int>>
 */
std::vector<std::vector<int>>
ExtractBezierPatchIDs(const std::vector<std::vector<int>>& knot_multiplicities,
                      const int* degrees,
                      const int* n_patches_per_para_dim);
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
      ExtractBezierPatchIDs(input.SplinepyKnotMultiplicities(),
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
