#pragma once

#include <array>
#include <type_traits>

// bezman
#include <bezman/src/rational_bezier_spline.hpp>
#include <bezman/src/point.hpp>

#include <splinepy/splines/helpers/scalar_type_wrapper.hpp>
#include <splinepy/splines/splinepy_base.hpp>

namespace splinepy::splines {

template<std::size_t para_dim, std::size_t dim>
using RationalBezierSpline = bezman::RationalBezierSpline<
    static_cast<std::size_t>(para_dim),
    std::conditional_t<(dim > 1), bezman::Point<dim>, double>,
    double>;

template<std::size_t para_dim, std::size_t dim>
class RationalBezier : public splinepy::splines::SplinepyBase,
                       public RationalBezierSpline<para_dim, dim> {
public:
  static constexpr int kParaDim = static_cast<int>(para_dim);
  static constexpr int kDim = static_cast<int>(para_dim);

  using SplinepyBase_ = typename splinepy::splines::SplinepyBase;
  using Base_ = RatioanlBezierSpline<para_dim, dim>;
  // alias to enable helper functions.
  using ParametricCoordinate_ = typename bezman::Point<para_dim, double>;
  using Coordinate_ = typename Base_::PhysicalPointType_;
  using Weight_ = typename Base_::ScalarType_;
  using Derivative_ = typename std::array<std::size_t, para_dim>;
  using Dimension_ = std::size_t;

  Base_ RawPtrInitHelper(const double* degrees, const double* control_points
                         const double* weights) {

    std::array<std::size_t, para_dim> bm_degrees{};
    std::size_t ncps{1};

    // formulate degrees
    for (std::size_t i{}; i < para_dim; ++i) {
      bm_degrees[i] = degrees[i];
      ncps *= degrees[i] + 1;
    }

    // formulate weighted control_points and weights.
    std::vector<Coordinate_> bm_weighted_control_points(ncps);
    std::vector<Weight_> bm_weights(ncps);
    for (std::size_t i{}; i < ncps; ++i) {
      // weights
      bm_weights[i] = weights[i];
      // weighted cps
      if constexpr (dim > 1) {
        for (std::size_t j = 0; j < dim; j++) {
          bm_weighted_control_points[i][j] = control_points[i * dim + j] * weights[i];
        }
      } else {
        bm_weighted_control_points[i] = control_points[i] * weights[i];
      }
    }
    return Base_(bm_degrees, bm_weighted_control_points, bm_weights);
  }

  // rawptr based ctor
  RationalBezier(const double* degrees, const double* control_points, const double* weights)
      : Base_(RawPtrInitHelper(degrees, control_points, weights)) {}
  // inherit ctor
  using Base_::Base_;

  // function wrapper, also for helper functions.
  constexpr auto operator()(const ParametricCoordinate_& query) const {
    return Base_::Evaluate(query);
  }

  constexpr auto operator()(const ParametricCoordinate_& query,
                            const Derivative_& order) const {
    return Base_::EvaluateDerivative(query, order);
  }

  constexpr auto ElevateDegree(const Dimension_ p_dim) {
    return Base_::OrderElevateAlongParametricDimension(p_dim);
  }

  // required implementations
  virtual constexpr int SplinepyParaDim() const { return kParaDim; }

  virtual constexpr int SplinepyDim() const { return kDim; }

  virtual std::string SplinepyWhatAmI() const {
    return "Bezier, parametric dimension: " + std::to_string(SplinepyParaDim())
           + ", physical dimension: " + std::to_string(SplinepyDim());
  }

  virtual int SplinepyNumberOfControlPoints() const {
    return static_cast<int>(Base_::control_points.size());
  }

  virtual void SplinepyCurrentProperties(
      double* degrees,
      std::vector<std::vector<double>>* knot_vectors /* untouched */,
      double* control_points,
      double* weights) const {

    // degrees
    for (std::size_t i{}; i < kParaDim; ++i) {
      degrees[i] = static_cast<double>(Base_::GetDegrees()[i]);
    }

    // control_points and weights
    const std::size_t ncps = Base_::control_points.size();
    for (std::size_t i{}; i < ncps; ++i) {
      const double w = Base_::GetWeights()[i];
      weights[i] = w;
      double inv_weight = static_cast<double>(1.) / w;
      if constexpr (dim > 1) {
        for (std::size_t j{}; j < kDim; ++j) {
          control_points[i * kDim + j] = Base_::GetWeightedControlPoints()[i][j] * inv_weight;
        }
      } else {
        control_points[i] = Base_::GetWeightedControlPoints()[i] * inv_weight;
      }
    }
  }

  virtual void SplinepyParametricBounds(double* para_bounds) const {
    for (std::size_t i{}; i < kParaDim; ++i) {
      // lower bounds
      para_bounds[i] = 0.;
      // upper bounds
      para_bounds[kParaDim + i] = 1.;
    }
  }

  virtual void SplinepyEvaluate(const double* para_coord,
                                double* evaluated) const {
    splinepy::splines::helpers::ScalarTypeEvaluate(*this,
                                                   para_coord,
                                                   evaluated);
  }
  virtual void SplinepyDerivative(const double* para_coord,
                                  const int* orders,
                                  double* derived) const {
    splinepy::splines::helpers::ScalarTypeDerivative(*this,
                                                     para_coord,
                                                     orders,
                                                     derived);
  }

  virtual void SplinepyElevateDegree(const int& p_dim) {
    splinepy::splines::helpers::ScalarTypeElevateDegree(*this, p_dim);
  }
}; /* class Bezier */

/// dynamic creation of templated BSpline
std::shared_ptr<SplinepyBase>
SplinepyBase::SplinepyCreateRationalBezier(const int para_dim,
                                   const int dim,
                                   const double* degrees,
                                   const double* control_points,
                                   const double* weights) {
  switch (para_dim) {
  case 1:
    switch (dim) {
    case 1:
      return std::make_shared<RationalBezier<1, 1>>(degrees, control_points, weights);
    case 2:
      return std::make_shared<RationalBezier<1, 2>>(degrees, control_points, weights);
    }
  case 2:
    switch (dim) {
    case 1:
      return std::make_shared<RationalBezier<2, 1>>(degrees, control_points, weights);
    case 2:
      return std::make_shared<RationalBezier<2, 2>>(degrees, control_points, weights);
    }
  }
}

} // namespace splinepy::splines
