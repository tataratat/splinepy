#pragma once

#include <array>

/* helper functions for all spline types */

namespace splinepy::splines {


/// Fills raw pointer with degrees. applicable for both nurbs and bspline
template<typename SplineType, typename DegreeType>
inline void FillDegrees(const SplineType& spline,
                        DegreeType* degree_ptr) {
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
template<typename SplineType, typename CountType>
inline CountType GetNumberOfControlPoints(const SplineType& spline) {
  return static_cast<CountType>(
      spline.GetVectorSpace().GetNumberOfCoordinates()
  );
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



} /* namespace splinepy::splines */
