#pragma once

#include <array>

#include "Sources/Utilities/index.hpp"

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

template<typename ResolutionType = int, typename SplineType>
inline std::array<ResolutionType, SplineType::kParaDim>
GetControlMeshResolutions(const SplineType& spline) {
  std::array<ResolutionType, SplineType::kParaDim> control_mesh_res;
  const auto& degrees = spline.GetDegrees();

  if constexpr (SplineType::kHasKnotVectors) {
    const auto& knot_vectors = spline.GetKnotVectors();
    for (int i{}; i < SplineType::kParaDim; ++i) {
      control_mesh_res[i] =
          static_cast<ResolutionType>(knot_vectors[i]->GetSize())
          - static_cast<ResolutionType>(degrees[i]) - 1;
    }
  } else {
    for (int i{}; i < SplineType::kParaDim; ++i) {
      control_mesh_res[i] = static_cast<ResolutionType>(degrees[i]) + 1;
    }
  }

  return control_mesh_res;
}

template<typename SplineType>
inline void GetGrevilleAbscissae(const SplineType& spline,
                                 double* greville_abscissae,
                                 const int& i_para_dim) {
  // Precompute values
  const auto& degrees = spline.GetDegrees();
  // Determine size using control point resolution
  const int cmr = [&]() {
    if constexpr (SplineType::kHasKnotVectors) {
      const auto& knot_vectors = spline.GetKnotVectors();
      return static_cast<int>(knot_vectors[i_para_dim]->GetSize())
             - static_cast<int>(degrees[i_para_dim]) - 1;

    } else {
      return static_cast<int>(degrees[i_para_dim]) + 1;
    }
  }();

  // Fill vector with points
  const double inv_factor = 1. / static_cast<double>(degrees[i_para_dim]);
  for (int j{}; j < cmr; ++j) {
    if constexpr (SplineType::kHasKnotVectors) {
      //
      using IndexType = typename splinelib::Index;
      const auto& knot_vectors = spline.GetKnotVectors()[i_para_dim];
      double factor{};
      for (int k{0}; k < degrees[i_para_dim]; ++k) {
        factor += knot_vectors->operator[](IndexType(k + j + 1));
      }
      greville_abscissae[j] = static_cast<double>(inv_factor * factor);
    } else {
      greville_abscissae[j] = static_cast<double>(inv_factor * j);
    }
  }
}

} // namespace splinepy::splines::helpers
