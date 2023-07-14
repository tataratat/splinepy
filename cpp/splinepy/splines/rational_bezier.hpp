#pragma once

#include <array>
#include <type_traits>

// bezman
#include <bezman/src/bezier_group.hpp>
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

template<std::size_t para_dim>
using RationalBezierSplineType =
    bezman::RationalBezierSpline<para_dim, bezman::VertexView<double>, double>;

/// @brief Rational Bezier (Spline)
/// @tparam para_dim
/// @tparam dim
template<std::size_t para_dim>
class RationalBezier : public SplinepyBase,
                       public RationalBezierSplineType<para_dim> {
public:
  /// @brief Parametric dimension
  static constexpr int kParaDim = static_cast<int>(para_dim);
  /// @brief True iff it is rational Bezier
  static constexpr bool kIsRational = true;
  /// @brief True iff it has knot vectors
  static constexpr bool kHasKnotVectors = false;

  using SplinepyBase_ = SplinepyBase;

  // self
  template<std::size_t s_para_dim>
  using SelfTemplate_ = RationalBezier<s_para_dim>;
  using SelfBoundary_ = RationalBezier<para_dim - 1>;

  // bezman
  using Base_ = RationalBezierSplineType<para_dim>;
  // alias to enable helper functions.
  using ParametricCoordinate_ = typename bezman::Point<para_dim, double>;
  using Degrees_ = typename std::array<std::size_t, para_dim>;
  // using Coordinate_ = typename Base_::PhysicalPointType_;
  //  keeping the notation similar to BSplineLib
  using VectorSpace_ = typename Base_::ControlPointsType_;
  using Coordinate_ = typename VectorSpace_::ContainerType_;
  using Coordinates_ = typename VectorSpace_::ContainerType_; // flat vector
  using Weight_ = typename Base_::ScalarType_;
  using Derivative_ = typename std::array<std::size_t, para_dim>;
  using Dimension_ = std::size_t;
  // advanced use
  using ProximityBase_ = proximity::ProximityBase;
  template<int dim>
  using Proximity_ =
      proximity::Proximity<RationalBezier<para_dim>, dim>;

  /// @brief Create base
  /// @param degrees
  /// @param control_points
  /// @param weights
  /// @return
  static Base_ CreateBase(const int* degrees,
                          const double* control_points,
                          const double* weights,
                          const int dim);

  /// @brief Construct a new Rational Bezier spline based on raw pointer
  /// @param degrees
  /// @param control_points
  /// @param weights
  RationalBezier(const int* degrees,
                 const double* control_points,
                 const double* weights,
                 const int dim)
      : Base_(CreateBase(degrees, control_points, weights, dim)) {}

  /// @brief  Base (copy) constructor
  /// @param rhs
  RationalBezier(const Base_& rhs) : Base_(rhs){};

  /// @brief Inherited constructor
  using Base_::Base_;

  /// @brief Function wrapper, also for helper functions. Evaluates query
  /// parametric coordinate
  /// @param query
  constexpr auto operator()(const ParametricCoordinate_& query) const;

  /// @brief Evaluates derivative of spline
  /// @param query
  /// @param order
  constexpr auto operator()(const ParametricCoordinate_& query,
                            const Derivative_& order) const;

  /// @brief Elevate degree
  /// @param p_dim
  constexpr auto ElevateDegree(const Dimension_ p_dim);

  // required implementations
  /// @brief Parametric dimension
  virtual int SplinepyParaDim() const { return kParaDim; }

  /// @brief Physical dimension
  virtual int SplinepyDim() const {
    return Base_::GetWeightedControlPoints().GetDimension();
  }

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
  virtual int SplinepyNumberOfSupports() const;

  /// @brief Current properties
  virtual void SplinepyCurrentProperties(
      int* degrees,
      std::vector<std::vector<double>>* knot_vectors /* untouched */,
      double* control_points,
      double* weights) const;

  virtual void SplinepyParametricBounds(double* para_bounds) const;

  virtual void SplinepyControlMeshResolutions(int* control_mesh_res) const;

  /// @brief Calculate Greville abscissae for Rational Bezier
  ///
  /// @param[out] greville_abscissae pointer to solution
  /// @param[in] i_para_dim parametric dimension
  virtual void SplinepyGrevilleAbscissae(double* greville_abscissae,
                                         const int& i_para_dim) const;

  virtual void SplinepyEvaluate(const double* para_coord,
                                double* evaluated) const;

  virtual void SplinepyDerivative(const double* para_coord,
                                  const int* orders,
                                  double* derived) const;

  virtual void SplinepyJacobian(const double* para_coord,
                                double* jacobians) const;

  virtual void SplinepyElevateDegree(const int& p_dim);

      virtual void SplinepyBasis(const double* para_coord, double* basis) const;

  virtual void SplinepyBasisDerivative(const double* para_coord,
                                       const int* order,
                                       double* basis_der) const;

  virtual void SplinepySupport(const double* para_coord, int* support) const;

  /// Basis Function values and their support IDs
  virtual void SplinepyBasisAndSupport(const double* para_coord,
                                       double* basis,
                                       int* support) const;

  /// Basis Function Derivative and their support IDs
  virtual void SplinepyBasisDerivativeAndSupport(const double* para_coord,
                                                 const int* orders,
                                                 double* basis_der,
                                                 int* support) const;

  virtual void SplinepyPlantNewKdTreeForProximity(const int* resolutions,
                                                  const int& nthreads);

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
                                        double* second_derivatives) const;

  /// only applicable to the splines of same para_dim, same type, and
  /// {1 or same} dim.
  virtual std::shared_ptr<SplinepyBase>
  SplinepyMultiply(const std::shared_ptr<SplinepyBase>& a) const;

  /// Spline addition.
  /// requires same type, para_dim, dim
  virtual std::shared_ptr<SplinepyBase>
  SplinepyAdd(const std::shared_ptr<SplinepyBase>& a) const;

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
      const std::shared_ptr<SplinepyBase>& inner_function) const;

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
      const std::shared_ptr<SplinepyBase>& inner_function) const;

  /// Spline composition.
  /// inner_function requirements:
  ///   1. Bezier Types
  ///   2. dim is same as outer_function's par_dim
  virtual std::shared_ptr<SplinepyBase>
  SplinepyCompose(const std::shared_ptr<SplinepyBase>& inner_function) const;

  virtual std::vector<std::shared_ptr<SplinepyBase>>
  SplinepySplit(const int& p_dim, const double& location) const;

  virtual std::shared_ptr<SplinepyBase>
  SplinepyDerivativeSpline(const int* orders) const;

  virtual std::shared_ptr<SplinepyBase>
  SplinepyExtractBoundary(const int& boundary_id);

  virtual std::vector<std::shared_ptr<SplinepyBase>>
  SplinepyExtractBezierPatches() const;

  virtual std::shared_ptr<SplinepyBase>
  SplinepyExtractDim(const int& phys_dim) const;

  /// Derivative of composition
  virtual std::shared_ptr<SplinepyBase> SplinepyCompositionDerivative(
      const std::shared_ptr<SplinepyBase>& inner,
      const std::shared_ptr<SplinepyBase>& inner_derivative) const;

  constexpr void CreateProximityIfNull();

  /// @brief Gets proximity
  constexpr ProximityBase_& GetProximity();

  /// @brief Gets proximity
  constexpr const ProximityBase_& GetProximity() const;

protected:
  /// @brief Shared pointer to proximity
  std::shared_ptr<ProximityBase_> proximity_ = nullptr;

}; /* class RationalBezier */

//#include <splinepy/splines/rational_bezier.inl>

} // namespace splinepy::splines

#include <splinepy/explicit/splinepy/rational_bezier_extern.hpp>
