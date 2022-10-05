#pragma once

#include <array>
#include <type_traits>

// bezman
#include <bezman/src/point.hpp>
#include <bezman/src/rational_bezier_spline.hpp>

#include <splinepy/splines/helpers/properties.hpp>
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
  static constexpr int kDim = static_cast<int>(dim);
  static constexpr bool kIsRational = true;
  static constexpr bool kHasKnotVectors = false;

  using SplinepyBase_ = typename splinepy::splines::SplinepyBase;
  using Base_ = RationalBezierSpline<para_dim, dim>;
  // alias to enable helper functions.
  using ParametricCoordinate_ = typename bezman::Point<para_dim, double>;
  using Coordinate_ = typename Base_::PhysicalPointType_;
  using Weight_ = typename Base_::ScalarType_;
  using Derivative_ = typename std::array<std::size_t, para_dim>;
  using Dimension_ = std::size_t;
  // advanced use
  using Proximity_ =
      splinepy::proximity::Proximity<RationalBezier<para_dim, dim>>;

  Base_ RawPtrInitHelper(const double* degrees,
                         const double* control_points,
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
          bm_weighted_control_points[i][j] =
              control_points[i * dim + j] * weights[i];
        }
      } else {
        bm_weighted_control_points[i] = control_points[i] * weights[i];
      }
    }
    return Base_(bm_degrees, bm_weighted_control_points, bm_weights);
  }

  // rawptr based ctor
  RationalBezier(const double* degrees,
                 const double* control_points,
                 const double* weights)
      : Base_(RawPtrInitHelper(degrees, control_points, weights)) {}
  // base (copy) ctor
  RationalBezier(const Base_& rhs) : Base_(rhs){};
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

  virtual std::string SplinepySplineName() const { return "RationalBezier"; }

  virtual std::string SplinepyWhatAmI() const {
    return "RationalBezier, parametric dimension: "
           + std::to_string(SplinepyParaDim())
           + ", physical dimension: " + std::to_string(SplinepyDim());
  }

  virtual bool SplinepyHasKnotVectors() const { return kHasKnotVectors; }

  virtual bool SplinepyIsRational() const { return kIsRational; }

  virtual int SplinepyNumberOfControlPoints() const {
    return static_cast<int>(Base_::GetWeightedControlPoints().size());
  }

  virtual int SplinepyNumberOfSupports() const {
    return splinepy::splines::helpers::GetNumberOfSupports(*this);
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
    const std::size_t ncps = Base_::GetWeightedControlPoints().size();
    for (std::size_t i{}; i < ncps; ++i) {
      const double w = Base_::GetWeights()[i];
      weights[i] = w;
      double inv_weight = static_cast<double>(1.) / w;
      if constexpr (dim > 1) {
        for (std::size_t j{}; j < kDim; ++j) {
          control_points[i * kDim + j] =
              Base_::GetWeightedControlPoints()[i][j] * inv_weight;
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

  /// only applicable to the splines of same para_dim, same type, and
  /// {1 or same} dim.
  virtual std::shared_ptr<SplinepyBase>
  SplinepyMultiply(const std::shared_ptr<SplinepyBase>& a) const {

    SplinepySplineNameMatches(
        *this,
        *a,
        "Spline multiplication requires splines of same type.",
        true);
    SplinepyParaDimMatches(
        *this,
        *a,
        "Spline multiplication requires splines of same parametric dimension.",
        true);
    // check dimension if a is not a scalar spline
    if (a->SplinepyDim() != 1) {
      SplinepyDimMatches(
          *this,
          *a,
          (std::string) "Spline multiplication requires splines of either 1 or "
              + "same physical dimension.",
          true);
    }

    // good to multiply
    if (a->SplinepyDim() == 1) {
      auto true_a = static_cast<RationalBezier<para_dim, 1>&>(*a);
      return std::make_shared<RationalBezier<para_dim, dim>>(
          this->Base_::operator*(true_a));
    } else {
      auto true_a = static_cast<RationalBezier<para_dim, dim>&>(*a);
      return std::make_shared<RationalBezier<para_dim, 1>>(
          this->Base_::operator*(true_a));
    }
  }

  /// Spline addition.
  /// requires same type, para_dim, dim
  virtual std::shared_ptr<SplinepyBase>
  SplinepyAdd(const std::shared_ptr<SplinepyBase>& a) const {
    SplinepySplineNameMatches(
        *this,
        *a,
        "Spline addition requires splines of the same time.",
        true);
    SplinepyParaDimMatches(
        *this,
        *a,
        "Spline addition requires splines of the same parametric dimension.",
        true);
    SplinepyDimMatches(
        *this,
        *a,
        "Spline addition requires splines of the same physical dimension",
        true);

    return std::make_shared<RationalBezier<para_dim, dim>>(
        this->Base_::operator+(
            static_cast<RationalBezier<para_dim, dim>&>(*a)));
  }

  /// Spline composition.
  /// inner_function requirements:
  ///   1. Bezier Types
  ///   2. dim is same as outer_function's par_dim
  virtual std::shared_ptr<SplinepyBase>
  SplinepyCompose(const std::shared_ptr<SplinepyBase>& inner_function) const {
    // type check
    if (inner_function->SplinepySplineName().find("Bezier")
        == std::string::npos) {
      splinepy::utils::PrintAndThrowError(
          "Bezier composition requires inner function to be a bezier type.",
          "Given inner function -",
          inner_function->SplinepyWhatAmI());
    }

    // composable?
    if (inner_function->SplinepyDim() != this->SplinepyParaDim()) {
      splinepy::utils::PrintAndThrowError(
          "Spline composition requires inner function to have same physical",
          "dimension as outer function's parametric dimension.",
          "Outer Function:",
          this->SplinepyWhatAmI(),
          "/",
          "Inner Function:",
          inner_function->SplinepyWhatAmI());
    }

    // compose - this time, it is always rational.
    if (inner_function->SplinepyIsRational()) {
      switch (inner_function->SplinepyParaDim()) {
      case 1:
        return std::make_shared<RationalBezier<1, dim>>(this->Compose(
            static_cast<RationalBezier<1, para_dim>&>(*inner_function)));
      case 2:
        return std::make_shared<RationalBezier<2, dim>>(this->Compose(
            static_cast<RationalBezier<2, para_dim>&>(*inner_function)));
      }
    } else {
      switch (inner_function->SplinepyParaDim()) {
      case 1:
        return std::make_shared<RationalBezier<1, dim>>(
            this->Compose(static_cast<Bezier<1, para_dim>&>(*inner_function)));
      case 2:
        return std::make_shared<RationalBezier<2, dim>>(
            this->Compose(static_cast<Bezier<2, para_dim>&>(*inner_function)));
      }
    }
  }

  virtual std::vector<std::shared_ptr<SplinepyBase>>
  SplinepySplit(const int& p_dim, const double& location) const {
    // split
    auto bm_splitted =
        Base_::SplitAtPosition(location, static_cast<std::size_t>(p_dim));

    // make it splinepybase
    std::vector<std::shared_ptr<SplinepyBase>> splitted;
    const std::size_t n_splitted = bm_splitted.size();
    splitted.reserve(n_splitted); // this should be always 2
    for (std::size_t i{}; i < n_splitted; ++i) {
      splitted.emplace_back(
          std::make_shared<RationalBezier<para_dim, dim>>(bm_splitted[i]));
    }

    return splitted;
  }

  virtual std::shared_ptr<SplinepyBase>
  SplinepyDerivativeSpline(const int* orders) const {
    // copy construct to start
    Base_ derived_bez{*this};
    // derive
    for (std::size_t i{}; i < para_dim; ++i) {
      for (int j{}; j < orders[i]; ++j) {
        derived_bez = derived_bez.DerivativeWRTParametricDimension(i);
      }
    }

    return std::make_shared<RationalBezier<para_dim, dim>>(derived_bez);
  }

  virtual std::vector<std::shared_ptr<SplinepyBase>>
  SplinepyExtractBezierPatches() const {
    // should copy
    return {std::make_shared<RationalBezier<para_dim, dim>>(*this)};
  }

  Proximity_& GetProximity() {
    if (!proximity_initialized_) {
      proximity_ = std::make_shared<Proximity_>(*this);
      proximity_initialized_ = true;
    }
    return *proximity_;
  }

protected:
  std::shared_ptr<Proximity_> proximity_;
  bool proximity_initialized_ = false;

}; /* class RationalBezier */

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
      return std::make_shared<RationalBezier<1, 1>>(degrees,
                                                    control_points,
                                                    weights);
    case 2:
      return std::make_shared<RationalBezier<1, 2>>(degrees,
                                                    control_points,
                                                    weights);
    }
  case 2:
    switch (dim) {
    case 1:
      return std::make_shared<RationalBezier<2, 1>>(degrees,
                                                    control_points,
                                                    weights);
    case 2:
      return std::make_shared<RationalBezier<2, 2>>(degrees,
                                                    control_points,
                                                    weights);
    }
  }
}

} // namespace splinepy::splines
