#pragma once

#include <array>
#include <type_traits>

// bezman
#include <bezman/src/point.hpp>
#include <bezman/src/rational_bezier_spline.hpp>

#include <splinepy/explicit/bezman/rational_bezier_extern.hpp>
#include <splinepy/proximity/proximity.hpp>
#include <splinepy/splines/bezier.hpp>
#include <splinepy/splines/helpers/basis_functions.hpp>
#include <splinepy/splines/helpers/extract.hpp>
#include <splinepy/splines/helpers/properties.hpp>
#include <splinepy/splines/helpers/scalar_type_wrapper.hpp>
#include <splinepy/splines/splinepy_base.hpp>

namespace splinepy::splines {

template<std::size_t para_dim, std::size_t dim>
using RationalBezierSplineType = bezman::RationalBezierSpline<
    static_cast<std::size_t>(para_dim),
    std::conditional_t<(dim > 1), bezman::Point<dim>, double>,
    double>;

/// @brief Rational Bezier (Spline)
/// @tparam para_dim
/// @tparam dim
template<std::size_t para_dim, std::size_t dim>
class RationalBezier : public splinepy::splines::SplinepyBase,
                       public RationalBezierSplineType<para_dim, dim> {
public:
  /// @brief Parametric dimension
  static constexpr int kParaDim = static_cast<int>(para_dim);
  /// @brief Physical dimension
  static constexpr int kDim = static_cast<int>(dim);
  /// @brief True iff it is rational Bezier
  static constexpr bool kIsRational = true;
  /// @brief True iff it has knot vectors
  static constexpr bool kHasKnotVectors = false;

  using SplinepyBase_ = splinepy::splines::SplinepyBase;
  using WeightedControlPointPointers_ =
      typename SplinepyBase_::WeightedControlPointPointers_;
  using WeightPointers_ = typename SplinepyBase_::WeightPointers_;

  // self
  template<std::size_t s_para_dim, std::size_t s_dim>
  using SelfTemplate_ = RationalBezier<s_para_dim, s_dim>;

  // bezman
  using Base_ = RationalBezierSplineType<para_dim, dim>;
  // alias to enable helper functions.
  using ParametricCoordinate_ = typename bezman::Point<para_dim, double>;
  using Degrees_ = typename std::array<std::size_t, para_dim>;
  using Coordinate_ = typename Base_::PhysicalPointType_;
  using Coordinates_ = typename std::vector<Coordinate_>;
  using Weight_ = typename Base_::ScalarType_;
  using Derivative_ = typename std::array<std::size_t, para_dim>;
  using Dimension_ = std::size_t;
  // advanced use
  using Proximity_ =
      splinepy::proximity::Proximity<RationalBezier<para_dim, dim>>;

  /// @brief Create base
  /// @param degrees
  /// @param control_points
  /// @param weights
  /// @return
  static Base_ CreateBase(const int* degrees,
                          const double* control_points,
                          const double* weights) {

    std::array<std::size_t, para_dim> bm_degrees{};
    std::size_t ncps{1};

    // formulate degrees
    for (std::size_t i{}; i < para_dim; ++i) {
      bm_degrees[i] = degrees[i];
      ncps *= degrees[i] + 1;
    }

    // formulate control_points and weights. bezman will apply weights
    std::vector<Coordinate_> bm_control_points(ncps);
    std::vector<Weight_> bm_weights(ncps);
    for (std::size_t i{}; i < ncps; ++i) {
      // weights
      bm_weights[i] = weights[i];
      // weighted cps
      if constexpr (dim > 1) {
        for (std::size_t j = 0; j < dim; j++) {
          bm_control_points[i][j] = control_points[i * dim + j];
        }
      } else {
        bm_control_points[i] = control_points[i];
      }
    }
    return Base_(bm_degrees, bm_control_points, bm_weights);
  }

  /// @brief Construct a new Rational Bezier spline based on raw pointer
  /// @param degrees
  /// @param control_points
  /// @param weights
  RationalBezier(const int* degrees,
                 const double* control_points,
                 const double* weights)
      : Base_(CreateBase(degrees, control_points, weights)) {}
  /// @brief  Base (copy) constructor
  /// @param rhs
  RationalBezier(const Base_& rhs) : Base_(rhs){};
  /// @brief Inherited constructor
  using Base_::Base_;

  /// @brief Function wrapper, also for helper functions. Evaluates query
  /// parametric coordinate
  /// @param query
  constexpr auto operator()(const ParametricCoordinate_& query) const {
    return Base_::Evaluate(query);
  }

  /// @brief Evaluates derivative of spline
  /// @param query
  /// @param order
  constexpr auto operator()(const ParametricCoordinate_& query,
                            const Derivative_& order) const {
    return Base_::EvaluateDerivative(query, order);
  }

  /// @brief Elevate degree
  /// @param p_dim
  constexpr auto ElevateDegree(const Dimension_ p_dim) {
    return Base_::OrderElevateAlongParametricDimension(p_dim);
  }

  // required implementations
  /// @brief Parametric dimension
  virtual int SplinepyParaDim() const { return kParaDim; }

  /// @brief Physical dimension
  virtual int SplinepyDim() const { return kDim; }

  /// @brief Name of Bezier curve
  virtual std::string SplinepySplineName() const { return "RationalBezier"; }

  /// @brief  What am I?
  /// @return String
  virtual std::string SplinepyWhatAmI() const {
    return "RationalBezier, parametric dimension: "
           + std::to_string(SplinepyParaDim())
           + ", physical dimension: " + std::to_string(SplinepyDim());
  }

  /// @brief Returns true iff Bezier curve has knot vectors
  virtual bool SplinepyHasKnotVectors() const { return kHasKnotVectors; }

  /// @brief Returns true iff Bezier curve is rational
  /// @return Bool
  virtual bool SplinepyIsRational() const { return kIsRational; }

  /// @brief Number of control points
  /// @return int
  virtual int SplinepyNumberOfControlPoints() const {
    return static_cast<int>(Base_::GetWeightedControlPoints().size());
  }

  /// @brief Number of supports
  /// @return int
  virtual int SplinepyNumberOfSupports() const {
    return splinepy::splines::helpers::GetNumberOfSupports(*this);
  }

  /// @brief Current properties
  virtual void SplinepyCurrentProperties(
      int* degrees,
      std::vector<std::vector<double>>* knot_vectors /* untouched */,
      double* control_points,
      double* weights) const {

    // degrees
    if (degrees) {
      for (std::size_t i{}; i < kParaDim; ++i) {
        degrees[i] = static_cast<int>(Base_::GetDegrees()[i]);
      }
    }

    // control_points and weights
    if (control_points || weights) {
      const std::size_t ncps = Base_::GetWeightedControlPoints().size();
      for (std::size_t i{}; i < ncps; ++i) {
        const double w = Base_::GetWeights()[i];
        if (weights) {
          weights[i] = w;
        }
        if (control_points) {
          double inv_weight = static_cast<double>(1.) / w;
          if constexpr (dim > 1) {
            for (std::size_t j{}; j < kDim; ++j) {
              control_points[i * kDim + j] =
                  Base_::GetWeightedControlPoints()[i][j] * inv_weight;
            }
          } else {
            control_points[i] =
                Base_::GetWeightedControlPoints()[i] * inv_weight;
          }
        }
      }
    }
  }

  virtual std::shared_ptr<WeightedControlPointPointers_>
  SplinepyWeightedControlPointPointers() {
    if (SplinepyBase_::control_point_pointers_) {
      return SplinepyBase_::control_point_pointers_;
    }

    // create weighted cps
    auto wcpp = std::make_shared<WeightedControlPointPointers_>();
    wcpp->dim_ = kDim;
    wcpp->for_rational_ = kIsRational;
    wcpp->coordinate_begins_.reserve(SplinepyNumberOfControlPoints());
    if constexpr (kDim == 1) {
      for (auto& w_control_point : Base_::GetWeightedControlPoints()) {
        wcpp->coordinate_begins_.push_back(&w_control_point);
      }
    } else {
      for (auto& w_control_point : Base_::GetWeightedControlPoints()) {
        wcpp->coordinate_begins_.push_back(w_control_point.data());
      }
    }

    // create weights
    auto w = std::make_shared<WeightPointers_>();
    w->weights_.reserve(SplinepyNumberOfControlPoints());
    for (auto& weight : Base_::GetWeights()) {
      w->weights_.push_back(&weight);
    }

    // reference each other
    w->control_point_pointers_ = wcpp; // weak_ptr
    wcpp->weight_pointers_ = w;

    SplinepyBase_::control_point_pointers_ = wcpp;

    return wcpp;
  }

  virtual std::shared_ptr<WeightPointers_> SplinepyWeightPointers() {
    return SplinepyWeightedControlPointPointers()->weight_pointers_;
  }

  virtual void SplinepyParametricBounds(double* para_bounds) const {
    for (std::size_t i{}; i < kParaDim; ++i) {
      // lower bounds
      para_bounds[i] = 0.;
      // upper bounds
      para_bounds[kParaDim + i] = 1.;
    }
  }

  virtual void SplinepyControlMeshResolutions(int* control_mesh_res) const {
    const auto cm_res =
        splinepy::splines::helpers::GetControlMeshResolutions(*this);
    std::copy_n(cm_res.begin(), para_dim, control_mesh_res);
  }

  /// @brief Calculate Greville abscissae for Rational Bezier
  ///
  /// @param[out] greville_abscissae pointer to solution
  /// @param[in] i_para_dim parametric dimension
  virtual void SplinepyGrevilleAbscissae(double* greville_abscissae,
                                         const int& i_para_dim) const {
    splinepy::splines::helpers::GetGrevilleAbscissae(*this,
                                                     greville_abscissae,
                                                     i_para_dim);
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

  virtual void SplinepyJacobian(const double* para_coord,
                                double* jacobians) const {
    splinepy::splines::helpers::ScalarTypeJacobian(*this,
                                                   para_coord,
                                                   jacobians);
  }

  virtual void SplinepyElevateDegree(const int& p_dim) {
    splinepy::splines::helpers::ScalarTypeElevateDegree(*this, p_dim);
  }

  virtual void SplinepyBasis(const double* para_coord, double* basis) const {
    splinepy::splines::helpers::BezierBasis(*this, para_coord, basis);
  }

  virtual void SplinepyBasisDerivative(const double* para_coord,
                                       const int* order,
                                       double* basis_der) const {
    splinepy::splines::helpers::BezierBasisDerivative(*this,
                                                      para_coord,
                                                      order,
                                                      basis_der);
  }

  virtual void SplinepySupport(const double* para_coord, int* support) const {
    splinepy::splines::helpers::BezierSupport(*this, para_coord, support);
  }

  /// Basis Function values and their support IDs
  virtual void SplinepyBasisAndSupport(const double* para_coord,
                                       double* basis,
                                       int* support) const {

    SplinepyBasis(para_coord, basis);
    SplinepySupport(para_coord, support);
  }

  /// Basis Function Derivative and their support IDs
  virtual void SplinepyBasisDerivativeAndSupport(const double* para_coord,
                                                 const int* orders,
                                                 double* basis_der,
                                                 int* support) const {
    SplinepyBasisDerivative(para_coord, orders, basis_der);
    SplinepySupport(para_coord, support);
  }

  virtual void SplinepyPlantNewKdTreeForProximity(const int* resolutions,
                                                  const int& nthreads) {
    splinepy::splines::helpers::ScalarTypePlantNewKdTreeForProximity(
        *this,
        resolutions,
        nthreads);
  }

  /// Verbose proximity query - make sure to plant a kdtree first.
  virtual void SplinepyVerboseProximity(const double* query,
                                        const double& tolerance,
                                        const int& max_iterations,
                                        const bool aggressive_bounds,
                                        double* para_coord,
                                        double* phys_coord,
                                        double* phys_diff,
                                        double& distance,
                                        double& convergence_norm,
                                        double* first_derivatives,
                                        double* second_derivatives) const {
    GetProximity().VerboseQuery(query,
                                tolerance,
                                max_iterations,
                                aggressive_bounds,
                                para_coord,
                                phys_coord,
                                phys_diff,
                                distance,
                                convergence_norm,
                                first_derivatives,
                                second_derivatives);
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

  /**
   * @brief Convert the composed splines into base splines for python side
   *
   * @tparam is_rational Flag for rational beziers
   * @tparam inner_para_dim parametric of inner function
   * @param inner_function inner function of composition
   * @return std::vector<std::shared_ptr<SplinepyBase>> Vector of SplinepyBase
   * splines
   */
  template<bool is_rational, size_t inner_para_dim>
  std::vector<std::shared_ptr<SplinepyBase>> ConvertComposeToBase(
      const std::shared_ptr<SplinepyBase>& inner_function) const {
    if constexpr (is_rational) {
      auto super_return = this->ComposeSensitivity(
          static_cast<RationalBezier<inner_para_dim, para_dim>&>(
              *inner_function));

      std::vector<std::shared_ptr<SplinepyBase>> to_return(super_return.size());

      for (size_t i{}; i < super_return.size(); ++i) {
        to_return[i] = std::make_shared<RationalBezier<inner_para_dim, 1>>(
            super_return[i]);
      }
      return to_return;
    } else {
      auto super_return = this->ComposeSensitivity(
          static_cast<Bezier<inner_para_dim, para_dim>&>(*inner_function));

      std::vector<std::shared_ptr<SplinepyBase>> to_return(super_return.size());

      for (size_t i{}; i < super_return.size(); ++i) {
        to_return[i] = std::make_shared<RationalBezier<inner_para_dim, 1>>(
            super_return[i]);
      }
      return to_return;
    }
  }

  /**
   * @brief Compute sensitivities of Composed splines with respect to outer
   * functions control points
   *
   *  Spline composition.
   *  inner_function requirements:
   *    1. Bezier Types
   *    2. dim is same as outer_function's par_dim
   *
   * @param inner_function inner function of composition
   * @return std::vector<std::shared_ptr<SplinepyBase>>
   */
  virtual std::vector<std::shared_ptr<SplinepyBase>>
  SplinepyComposeSensitivities(
      const std::shared_ptr<SplinepyBase>& inner_function) const {
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

    // compose - return correct type.
    if (inner_function->SplinepyIsRational()) {
      switch (inner_function->SplinepyParaDim()) {
      case 1:
        return ConvertComposeToBase<true, 1>(inner_function);
      case 2:
        return ConvertComposeToBase<true, 2>(inner_function);
      case 3:
        return ConvertComposeToBase<true, 3>(inner_function);
#ifdef SPLINEPY_MORE
      case 4:
        return ConvertComposeToBase<true, 4>(inner_function);
      case 5:
        return ConvertComposeToBase<true, 5>(inner_function);
      case 6:
        return ConvertComposeToBase<true, 6>(inner_function);
      case 7:
        return ConvertComposeToBase<true, 7>(inner_function);
      case 8:
        return ConvertComposeToBase<true, 8>(inner_function);
      case 9:
        return ConvertComposeToBase<true, 9>(inner_function);
      case 10:
        return ConvertComposeToBase<true, 10>(inner_function);
#endif
      default:
        splinepy::utils::PrintAndThrowError(
            "Something went wrong during ComposeSensitivity. Please help us by "
            "writing an issue about this case at [ "
            "github.com/tataratat/splinepy ]");
      }
    } else {
      switch (inner_function->SplinepyParaDim()) {
      case 1:
        return ConvertComposeToBase<false, 1>(inner_function);
      case 2:
        return ConvertComposeToBase<false, 2>(inner_function);
      case 3:
        return ConvertComposeToBase<false, 3>(inner_function);
#ifdef SPLINEPY_MORE
      case 4:
        return ConvertComposeToBase<false, 4>(inner_function);
      case 5:
        return ConvertComposeToBase<false, 5>(inner_function);
      case 6:
        return ConvertComposeToBase<false, 6>(inner_function);
      case 7:
        return ConvertComposeToBase<false, 7>(inner_function);
      case 8:
        return ConvertComposeToBase<false, 8>(inner_function);
      case 9:
        return ConvertComposeToBase<false, 9>(inner_function);
      case 10:
        return ConvertComposeToBase<false, 10>(inner_function);
#endif
      default:
        splinepy::utils::PrintAndThrowError(
            "Something went wrong during ComposeSensitivity. Please help us by "
            "writing an issue about this case at [ "
            "github.com/tataratat/splinepy ]");
      }
    }
    splinepy::utils::PrintAndThrowError(
        "Something went wrong during ComposeSensitivity. Please help us by "
        "writing an issue about this case at [ "
        "github.com/tataratat/splinepy ]");
    // make compiler happy
    return std::vector<std::shared_ptr<SplinepyBase>>{};
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
      case 3:
        return std::make_shared<RationalBezier<3, dim>>(this->Compose(
            static_cast<RationalBezier<3, para_dim>&>(*inner_function)));
#ifdef SPLINEPY_MORE
      case 4:
        return std::make_shared<RationalBezier<4, dim>>(this->Compose(
            static_cast<RationalBezier<4, para_dim>&>(*inner_function)));
      case 5:
        return std::make_shared<RationalBezier<5, dim>>(this->Compose(
            static_cast<RationalBezier<5, para_dim>&>(*inner_function)));
      case 6:
        return std::make_shared<RationalBezier<6, dim>>(this->Compose(
            static_cast<RationalBezier<6, para_dim>&>(*inner_function)));
      case 7:
        return std::make_shared<RationalBezier<7, dim>>(this->Compose(
            static_cast<RationalBezier<7, para_dim>&>(*inner_function)));
      case 8:
        return std::make_shared<RationalBezier<8, dim>>(this->Compose(
            static_cast<RationalBezier<8, para_dim>&>(*inner_function)));
      case 9:
        return std::make_shared<RationalBezier<9, dim>>(this->Compose(
            static_cast<RationalBezier<9, para_dim>&>(*inner_function)));
      case 10:
        return std::make_shared<RationalBezier<10, dim>>(this->Compose(
            static_cast<RationalBezier<10, para_dim>&>(*inner_function)));
#endif
      default:
        splinepy::utils::PrintAndThrowError(
            "Something went wrong during Compose. Please help us by writing an "
            "issue about this case at [ github.com/tataratat/splinepy ]");
      }
    } else {
      switch (inner_function->SplinepyParaDim()) {
      case 1:
        return std::make_shared<RationalBezier<1, dim>>(
            this->Compose(static_cast<Bezier<1, para_dim>&>(*inner_function)));
      case 2:
        return std::make_shared<RationalBezier<2, dim>>(
            this->Compose(static_cast<Bezier<2, para_dim>&>(*inner_function)));
      case 3:
        return std::make_shared<RationalBezier<3, dim>>(
            this->Compose(static_cast<Bezier<3, para_dim>&>(*inner_function)));
#ifdef SPLINEPY_MORE
      case 4:
        return std::make_shared<RationalBezier<4, dim>>(
            this->Compose(static_cast<Bezier<4, para_dim>&>(*inner_function)));
      case 5:
        return std::make_shared<RationalBezier<5, dim>>(
            this->Compose(static_cast<Bezier<5, para_dim>&>(*inner_function)));
      case 6:
        return std::make_shared<RationalBezier<6, dim>>(
            this->Compose(static_cast<Bezier<6, para_dim>&>(*inner_function)));
      case 7:
        return std::make_shared<RationalBezier<7, dim>>(
            this->Compose(static_cast<Bezier<7, para_dim>&>(*inner_function)));
      case 8:
        return std::make_shared<RationalBezier<8, dim>>(
            this->Compose(static_cast<Bezier<8, para_dim>&>(*inner_function)));
      case 9:
        return std::make_shared<RationalBezier<9, dim>>(
            this->Compose(static_cast<Bezier<9, para_dim>&>(*inner_function)));
      case 10:
        return std::make_shared<RationalBezier<10, dim>>(
            this->Compose(static_cast<Bezier<10, para_dim>&>(*inner_function)));
#endif
      default:
        splinepy::utils::PrintAndThrowError(
            "Something went wrong during Compose. Please help us by writing an "
            "issue about this case at [ github.com/tataratat/splinepy ]");
      }
    }
    splinepy::utils::PrintAndThrowError(
        "Something went very wrong during Compose. Please help us by writing"
        "an issue about this case at [ github.com/tataratat/splinepy ]");
    // make compiler happy
    return std::shared_ptr<SplinepyBase>{};
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

  virtual std::shared_ptr<SplinepyBase>
  SplinepyExtractBoundary(const int& boundary_id) {
    return splinepy::splines::helpers::ExtractBoundaryMeshSlice(*this,
                                                                boundary_id);
  }

  virtual std::vector<std::shared_ptr<SplinepyBase>>
  SplinepyExtractBezierPatches() const {
    // should copy
    return {std::make_shared<RationalBezier<para_dim, dim>>(*this)};
  }

  virtual std::shared_ptr<SplinepyBase>
  SplinepyExtractDim(const int& phys_dim) const {
    return std::make_shared<RationalBezier<para_dim, 1>>(
        Base_::ExtractDimension(static_cast<std::size_t>(phys_dim)));
  }

  /// Derivative of composition
  virtual std::shared_ptr<SplinepyBase> SplinepyCompositionDerivative(
      const std::shared_ptr<SplinepyBase>& inner,
      const std::shared_ptr<SplinepyBase>& inner_derivative) const {
    // type check
    if (inner->SplinepySplineName().find("Bezier") == std::string::npos) {
      splinepy::utils::PrintAndThrowError(
          "Bezier composition derivative requires inner function to be",
          "a bezier type. Given inner function -",
          inner->SplinepyWhatAmI());
    }

    // composable?
    if (inner->SplinepyDim() != this->SplinepyParaDim()) {
      splinepy::utils::PrintAndThrowError(
          "Spline composition requires inner function to have same physical",
          "dimension as outer function's parametric dimension.",
          "Outer Function:",
          this->SplinepyWhatAmI(),
          "/",
          "Inner Function:",
          inner->SplinepyWhatAmI());
    }
    std::array<int, para_dim> order_query;
    order_query.fill(0);
    order_query[0] = 1;
    std::shared_ptr<SplinepyBase> composition_derivative =
        this->SplinepyDerivativeSpline(order_query.data())
            ->SplinepyCompose(inner)
            ->SplinepyMultiply(inner_derivative->SplinepyExtractDim(0));

    for (std::size_t i{1}; i < para_dim; ++i) {
      order_query.fill(0);
      order_query[i] = 1;
      composition_derivative = composition_derivative->SplinepyAdd(
          this->SplinepyDerivativeSpline(order_query.data())
              ->SplinepyCompose(inner)
              ->SplinepyMultiply(inner_derivative->SplinepyExtractDim(i)));
    }
    return composition_derivative;
  }

  /// @brief Get proximity
  constexpr Proximity_& GetProximity() { return *proximity_; }
  /// @brief Get proximity
  constexpr const Proximity_& GetProximity() const { return *proximity_; }

protected:
  /// @brief Shared pointer to proximity
  std::shared_ptr<Proximity_> proximity_ = std::make_shared<Proximity_>(*this);

}; /* class RationalBezier */

} // namespace splinepy::splines

#include <splinepy/explicit/splinepy/rational_bezier_extern.hpp>
