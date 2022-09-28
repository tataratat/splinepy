#pragma once

#include <array>

namespace splinepy::splines::helpers {

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

} // namespace splinepy::splines::helpers
