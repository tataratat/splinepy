#pragma once

namespace splinepy::splines {

// Forward declarations
template<int para_dim, int dim>
class Nurbs;

template<int para_dim, int dim>
class BSpline;

// Type Identifier
template<typename SplineType>
struct isNurbs {
  constexpr static bool value = false;
};

template<int para_dim, int dim>
struct isNurbs<Nurbs<para_dim, dim>> {
  constexpr static bool value = true;
};

template<typename SplineType>
constexpr bool isNurbs_v = isNurbs<SplineType>::value;

template<typename SplineType>
struct isBSpline {
  constexpr static bool value = false;
};

template<int para_dim, int dim>
struct isBSpline<BSpline<para_dim, dim>> {
  constexpr static bool value = true;
};

template<typename SplineType>
constexpr bool isBSpline_v = isBSpline<SplineType>::value;

} // namespace splinepy::splines
