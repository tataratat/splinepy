#pragma once

#include <splinepy/splines/splinepy_base.hpp>

namespace splinepy::splines::helpers {


/// returns boundary spline, which has one less para_dim.
template<typename SplineType>
std::shared_ptr<splinepy::splines::SplinepyBase>
ExtractBoundarySpline(SplineType& spline,
                      const int& plane_normal_axis,
                      const int& extrema) {
  if constexpr (SplineType::kParaDim == 1) {
    return std::shared_ptr<
    typename SplineType::template SelfTemplate_<SplineType::kParaDim - 1,
                                                SplineType::kDim>>{};

  } else {

  auto& grid_slicer = spline.GetControlMeshSampler();
  const auto ids_on_boundary =
      grid_slicer.GridPointIdsOnBoundary(plane_normal_axis, extrema);
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
    using KnotVector = typename KnotVectors::element_type;
    using VSpace = typename SelfBoundary::VectorSpace_;
    using Coords = typename VSpace::Coordinates_;

    // PSpace - let it create a new basis. Otherwise, it will make
    // pure copies of all basis function, which is not high degree friendly.
    Degrees b_degrees;
    KnotVectors b_knot_vectors;
    for (int i{}; i < SplineType::kParaDim; ++i) {
      if (i == plane_normal_axis) {
        continue;
      }
      b_degrees[i] = spline.GetDegrees()[i];
      // knotvectors are shared_ptrs. make sure this copies
      b_knot_vectors[i] =
          std::make_shared<KnotVector>(*spline.GetKnotVectors()[i]);
    }
    // create parametric space
    auto pspace = std::make_shared<PSpace>(b_degrees, b_knot_vectors);

    // VSpace - maybe it is also
    // create (weighted) vector space
    std::shared_ptr<VSpace> vspace{};
    auto& coords = vspace->GetCoordinates();
    coords.reserve(ids_on_boundary.size());
    const auto& rhs_coordinates = spline.GetCoordinates();
    for (const auto& id : ids_on_boundary) {
      if constexpr (SplineType::kIsRational) {
        // these are homogeneous coords
        coords.push_back(rhs_coordinates[id]);
      } else {
        coords.push_back(rhs_coordinates[id]);
      }
    }
    // assign boundary spline - uses base' ctor
    boundary_spline = std::make_shared<SelfBoundary>(pspace, vspace);
  } else {
    // hola bezier families
    // form degrees
    using Degrees = typename SelfBoundary::Degrees_;
    Degrees b_degrees;
    std::size_t ncps{1};
    for (int i{}; i < SplineType::kParaDim; ++i) {
      if (i == plane_normal_axis) {
        continue;
      }
      const auto& ds = spline.GetDegrees();
      const auto& d = ds[i];
      b_degrees[i] = d;
      ncps *= static_cast<std::size_t>(d + 1);
    }

    if constexpr (SplineType::kIsRational) {
      // only way to avoid unweight -> reweight is to init first with degrees
      // assign boundary spline
      boundary_spline = std::make_shared<SelfBoundary>(b_degrees);
      auto& b_cps = boundary_spline->GetWeightedControlPoints();
      auto& b_ws = boundary_spline->GetWeights();

      // copy
      for (std::size_t i{}; i < ncps; ++i) {
        const auto& id = ids_on_boundary[i];
        b_cps[i] = spline.GetWeightedControlPoints()[id];
        b_ws[i] = spline.GetWeights()[id];
      }
    } else {
      // non-rational bezier
      using Coordinates = typename SelfBoundary::Coordinates_;
      Coordinates b_coords;
      b_coords.reserve(ncps);
      for (const auto& id : ids_on_boundary) {
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

} /* namespace splinepy::splines::helpers */
