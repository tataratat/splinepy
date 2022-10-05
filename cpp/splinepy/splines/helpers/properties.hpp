#pragma once

#include <array>

namespace splinepy::splines::helpers {

/// Returns parametric bounds as (para_dim x 2) array.
template<typename SplineType>
inline std::array<std::array<double, SplineType::kParaDim>, 2>
GetParametricBounds(const SplineType& spline) {
  std::array<std::array<double, SplineType::kParaDim>, 2> parametric_bounds{};

  if constexpr (!SplineType::kHasKnotVectors) {
    parametric_bounds[0].fill(0.);
    parametric_bounds[1].fill(1.);
  } else {
    const auto& parameter_space = spline.GetParameterSpace();
    size_t i{};
    for (const auto& knotvector : parameter_space.GetKnotVectors()) {
      parametric_bounds[0][i] = knotvector->GetFront().Get();
      parametric_bounds[1][i] = knotvector->GetBack().Get();
      ++i;
    }
  }

  return parametric_bounds;
}

template<typename SplineType>
inline int GetNumberOfSupports(const SplineType& spline) {
  const auto& degrees = spline.GetDegrees();
  int n_supports{1};
  for (int i{}; i < SplineType::kParaDim; ++i) {
    n_supports *= static_cast<int>(degrees[i]) + 1;
  }
  return n_supports;
}

} // namespace splinepy::splines::helpers
