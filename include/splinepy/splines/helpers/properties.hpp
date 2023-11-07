#pragma once

#include <array>

#include "BSplineLib/Utilities/index.hpp"

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
    const auto& degrees = parameter_space.GetDegrees();

    size_t i{};
    for (const auto& knotvector : parameter_space.GetKnotVectors()) {
      const int& degree_value = degrees[i];
      parametric_bounds[0][i] = knotvector->operator[](degree_value);
      parametric_bounds[1][i] =
          knotvector->operator[](knotvector->GetSize() - degree_value - 1);
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

/// @brief Computes Greville Abscissae
/// @tparam SplineType
/// @param[in] spline Input Splines
/// @param[out] greville_abscissae Output 1D array
/// @param[in] i_para_dim parametric dimension along which greville abscissaes
///                       are computed
/// @param[in] duplicate_tolerance if negative two greville abscissae can be
///                                equal, positive tolerance to avoid
///                                duplication of greville abscissae. Made to
///                                comply with C^(-1) splines. Tolerance
///                                represents difference between two greville
///                                abscissae for them to be considered equal
template<typename SplineType>
inline void GetGrevilleAbscissae(const SplineType& spline,
                                 double* greville_abscissae,
                                 const int& i_para_dim,
                                 const double& duplicate_tolerance) {
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
      const auto& knot_vectors = spline.GetKnotVectors()[i_para_dim];
      double factor{};
      for (int k{0}; k < degrees[i_para_dim]; ++k) {
        factor +=
            knot_vectors->operator[](typename bsplinelib::Index(k + j + 1));
      }
      greville_abscissae[j] = static_cast<double>(inv_factor * factor);
    } else {
      greville_abscissae[j] = static_cast<double>(inv_factor * j);
    }
  }

  // There can be no duplicates for bezier types, but if C^(-1), duplicates in
  // the knot vector result in duplicate greville abscissae which can lead to
  // problems in further computations, solved by mean filtering with neighbors
  if constexpr (SplineType::kHasKnotVectors) {
    if (duplicate_tolerance > 0.0) {
      double previous_knot{greville_abscissae[1]};
      for (int j{2}; j < cmr - 1; ++j) {
        if (std::abs(previous_knot - greville_abscissae[j])
            < duplicate_tolerance) {
          // @todo make a dynamic tolerance
          greville_abscissae[j - 1] =
              0.5 * (greville_abscissae[j - 2] + greville_abscissae[j - 1]);
          greville_abscissae[j] =
              0.5 * (greville_abscissae[j] + greville_abscissae[j + 1]);
        } else {
          previous_knot = greville_abscissae[j];
        }
      }
    }
  }
}

} // namespace splinepy::splines::helpers
