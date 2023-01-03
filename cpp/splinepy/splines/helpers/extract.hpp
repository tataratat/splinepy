#pragma once

#include <splinepy/splines/helpers/properties.hpp>
#include <splinepy/splines/splinepy_base.hpp>
#include <splinepy/utils/grid_points.hpp>
#include <splinepy/utils/print.hpp>

namespace splinepy::splines::helpers {

/// Determine boundary ID from axis and extrem-value
inline int ExtractBoundaryFromAxisAndExtrema(const int& axis,
                                             const int& extreme) {
  // Determine corresponding ID
  return (extreme > 0) ? 2 * axis + 1 : 2 * axis;
}

template<typename SplineType, typename IndexT>
std::shared_ptr<splinepy::splines::SplinepyBase>
ExtractControlMeshSliceFromIDs(const SplineType& spline,
                               const std::vector<IndexT>& indices,
                               const int& plane_normal_axis) {
  // Perform sanity check
  if constexpr (SplineType::kParaDim == 1) {
    splinepy::utils::PrintWarning(
        "Sorry, we don't support control mesh slicing"
        "of 1-Parametric Dim splines. Returning empty spline.");
    return std::shared_ptr<
        typename SplineType::template SelfTemplate_<SplineType::kParaDim,
                                                    SplineType::kDim>>{};
  } else {
    // boundary spline type
    using SelfBoundary =
        typename SplineType::template SelfTemplate_<SplineType::kParaDim - 1,
                                                    SplineType::kDim>;
    // prepare return
    std::shared_ptr<SelfBoundary> boundary_spline;

    // create new one
    if constexpr (SplineType::kHasKnotVectors) {
      // hola BSpline families
      using PSpace = typename SelfBoundary::ParameterSpace_;
      using Degrees = typename PSpace::Degrees_;
      using KnotVectors = typename PSpace::KnotVectors_;
      using KnotVector = typename KnotVectors::value_type::element_type;
      using VSpace = typename SelfBoundary::PhysicalSpace_;
      // PSpace - let it create a new basis. Otherwise, it will make
      // pure copies of all basis function, which is not high degree friendly.
      Degrees b_degrees;
      KnotVectors b_knot_vectors;
      int counter{};
      for (int i{}; i < SplineType::kParaDim; ++i) {
        if (i == plane_normal_axis) {
          continue;
        }
        b_degrees[counter] = spline.GetDegrees()[i];
        // knotvectors are shared_ptrs. make sure this copies
        b_knot_vectors[counter] =
            std::make_shared<KnotVector>(*spline.GetKnotVectors()[i]);
        counter++;
      }
      // create parametric space
      auto pspace = std::make_shared<PSpace>(b_knot_vectors, b_degrees);

      // VSpace - maybe it is also
      // create (weighted) vector space
      auto vspace = std::make_shared<VSpace>();
      auto& coords = vspace->GetCoordinates();
      coords.reserve(indices.size());
      const auto& rhs_coordinates = spline.GetCoordinates();
      for (const auto& id : indices) {
        coords.push_back(rhs_coordinates[id]);
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

template<bool switch_plane_id_to_extrema = false, typename SplineType>
std::shared_ptr<splinepy::splines::SplinepyBase>
ExtractBoundaryMeshSlice(const SplineType& spline, const int& boundary_id) {
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
        splinepy::splines::helpers::template GetControlMeshResolutions<
            std::size_t>(spline);
    const int plane_normal_axis = static_cast<int>(boundary_id / 2);
    const int plane_id =
        ((boundary_id % 2) == 0) ? 0 : cmr[plane_normal_axis] - 1;
    const auto ids_on_boundary = splinepy::utils::GridPoints<
        std::size_t,
        std::size_t,
        SplineType::kParaDim>::GridPointIdsOnHyperPlane(cmr,
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
        splinepy::splines::helpers::template GetControlMeshResolutions<
            std::size_t>(spline);
    const auto ids_on_boundary = splinepy::utils::GridPoints<
        std::size_t,
        std::size_t,
        SplineType::kParaDim>::GridPointIdsOnHyperPlane(cmr,
                                                        plane_normal_axis,
                                                        plane_id);
    return ExtractControlMeshSliceFromIDs(spline,
                                          ids_on_boundary,
                                          plane_normal_axis);
  }
}

} /* namespace splinepy::splines::helpers */
