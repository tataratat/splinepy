#pragma once

namespace splinepy::splines {

// Forward declarations
template<int para_dim, int dim>
class Nurbs;

template<int para_dim, int dim>
class BSpline;

// Type Identifier
/// @struct Nurbs identifier
/// @brief Nurbs identifier
/// @tparam SplineType
template<typename SplineType>
struct isNurbs {
  /// @brief Is true iff spline is Nurbs
  constexpr static bool value = false;
};

/// @struct Nurbs identifier
/// @brief Nurbs identifier
/// @tparam para_dim Parametric dimension
/// @tparam dim Physical dimension
template<int para_dim, int dim>
struct isNurbs<Nurbs<para_dim, dim>> {
  /// @brief Is true iff spline is Nurbs
  constexpr static bool value = true;
};

template<typename SplineType>
constexpr bool isNurbs_v = isNurbs<SplineType>::value;

/// @struct BSpline identifier
/// @brief BSpline identifier
/// @tparam SplineType
template<typename SplineType>
struct isBSpline {
  /// @brief Is true iff spline is BSpline
  constexpr static bool value = false;
};

/// @struct BSpline identifier
/// @brief BSpline identifier
/// @tparam para_dim Parametric dimension
/// @tparam dim Physical dimension
template<int para_dim, int dim>
struct isBSpline<BSpline<para_dim, dim>> {
  /// @brief Is true iff spline is BSpline
  constexpr static bool value = true;
};

template<typename SplineType>
constexpr bool isBSpline_v = isBSpline<SplineType>::value;

} // namespace splinepy::splines
