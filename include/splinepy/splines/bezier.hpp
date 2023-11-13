#pragma once

#include <array>
#include <type_traits>

// bezman
#include <bezman/src/bezier_spline.hpp>
#include <bezman/src/point.hpp>

#include "splinepy/proximity/proximity.hpp"
#include "splinepy/splines/splinepy_base.hpp"

namespace splinepy::splines {

template<std::size_t para_dim, std::size_t dim>
using BezierSplineType = bezman::BezierSpline<
    static_cast<std::size_t>(para_dim),
    std::conditional_t<(dim > 1), bezman::Point<dim>, double>,
    double>;

/// @brief Bezier class
/// @tparam para_dim Dimension of parametric space
/// @tparam dim Dimension of physical space
template<std::size_t para_dim, std::size_t dim>
class Bezier : public splinepy::splines::SplinepyBase,
               public BezierSplineType<para_dim, dim> {
public:
  /// @brief Parametric space dimension
  static constexpr int kParaDim = static_cast<int>(para_dim);
  /// @brief Physical space dimension
  static constexpr int kDim = static_cast<int>(dim);
  /// @brief True iff it is rational Bezier
  static constexpr bool kIsRational = false;
  /// @brief True iff spline has knot vectors
  static constexpr bool kHasKnotVectors = false;

  using SplinepyBase_ = splinepy::splines::SplinepyBase;
  template<std::size_t s_para_dim, std::size_t s_dim>
  using SelfTemplate_ = Bezier<s_para_dim, s_dim>;
  using BoundaryType_ = Bezier<para_dim - 1, dim>;
  using Base_ = BezierSplineType<para_dim, dim>;
  template<std::size_t b_para_dim, std::size_t b_dim>
  using BaseTemplate_ = BezierSplineType<b_para_dim, b_dim>;
  using ControlPointPointers_ = typename SplinepyBase_::ControlPointPointers_;
  // alias to enable helper functions.
  using ParametricCoordinate_ = typename bezman::Point<para_dim, double>;
  using Degrees_ = typename std::array<std::size_t, para_dim>;
  using Coordinate_ = typename Base_::PhysicalPointType_;
  using Coordinates_ = typename std::vector<Coordinate_>;
  using Derivative_ = typename std::array<std::size_t, para_dim>;
  using Dimension_ = std::size_t;
  using Proximity_ = splinepy::proximity::Proximity;

  /// @brief Creates Base for Bezier
  /// @param degrees
  /// @param control_points
  static Base_ CreateBase(const int* degrees, const double* control_points);

  /// @brief Constructor based on raw pointer
  /// @param degrees Degrees of the Bezier
  /// @param control_points Control points of the Bezier
  Bezier(const int* degrees, const double* control_points)
      : Base_(CreateBase(degrees, control_points)) {}
  /// @brief (Copy) constructor
  Bezier(const Base_& rhs) : Base_(rhs) {}
  /// @brief Inherited constructor
  using Base_::Base_;

  // function wrapper, also for helper functions.
  /// @brief Evaluates Bezier
  /// @param query
  constexpr auto operator()(const ParametricCoordinate_& query) const;
  /// @brief Evaluates derivative
  /// @param query
  /// @param order
  constexpr auto operator()(const ParametricCoordinate_& query,
                            const Derivative_& order) const;

  /// @brief Elevates degrees
  /// @param p_dim
  constexpr auto ElevateDegree(const Dimension_ p_dim);

  // required implementations
  /// @copydoc splinepy::splines::SplinepyBase::SplinepyParaDim
  virtual int SplinepyParaDim() const { return kParaDim; }

  /// @copydoc splinepy::splines::SplinepyBase::SplinepyDim
  virtual int SplinepyDim() const { return kDim; }

  /// @copydoc splinepy::splines::SplinepyBase::SplinepySplineName
  virtual std::string SplinepySplineName() const { return "Bezier"; }

  /// @copydoc splinepy::splines::SplinepyBase::SplinepyWhatAmI
  virtual std::string SplinepyWhatAmI() const {
    return "Bezier, parametric dimension: " + std::to_string(SplinepyParaDim())
           + ", physical dimension: " + std::to_string(SplinepyDim());
  }

  /// @copydoc splinepy::splines::SplinepyBase::SplinepyHasKnotVectors
  virtual bool SplinepyHasKnotVectors() const { return kHasKnotVectors; }

  /// @copydoc splinepy::splines::SplinepyBase::SplinepyIsRational
  virtual bool SplinepyIsRational() const { return kIsRational; }

  /// @copydoc splinepy::splines::SplinepyBase::SplinepyNumberOfControlPoints
  virtual int SplinepyNumberOfControlPoints() const;

  /// @copydoc splinepy::splines::SplinepyBase::SplinepyNumberOfSupports
  virtual int SplinepyNumberOfSupports() const;

  /// @copydoc splinepy::splines::SplinepyBase::SplinepyCurrentProperties
  virtual void SplinepyCurrentProperties(
      int* degrees,
      std::vector<std::vector<double>>* knot_vectors /* untouched */,
      double* control_points,
      double* weights /* untouched */) const;

  virtual std::shared_ptr<ControlPointPointers_> SplinepyControlPointPointers();

  virtual void SplinepyParametricBounds(double* para_bounds) const;

  virtual void SplinepyControlMeshResolutions(int* control_mesh_res) const;

  /// @brief Calculate Greville abscissae for Bezier
  ///
  /// @param[out] greville_abscissae pointer to solution
  /// @param[in] i_para_dim parametric dimension
  /// @param[in] duplicate_tolerance if negative two greville abscissae can be
  ///                                equal, positive tolerance to avoid
  ///                                duplication of greville abscissae. Made to
  ///                                comply with C^(-1) splines. Tolerance
  ///                                represents difference between two greville
  ///                                abscissae for them to be considered equal
  virtual void
  SplinepyGrevilleAbscissae(double* greville_abscissae,
                            const int& i_para_dim,
                            const double& duplicate_tolerance) const;

  virtual void SplinepyEvaluate(const double* para_coord,
                                double* evaluated) const;

  virtual void SplinepyDerivative(const double* para_coord,
                                  const int* orders,
                                  double* derived) const;
  virtual void SplinepyJacobian(const double* para_coord,
                                double* jacobians) const;
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

  /// only applicable to the splines of same para_dim, same type
  /// and {1, same} physical dim.
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

  virtual std::vector<std::shared_ptr<SplinepyBase>>
  SplinepyExtractBezierPatches() const;

  virtual std::shared_ptr<SplinepyBase>
  SplinepyExtractBoundary(const int& boundary_id);

  virtual std::shared_ptr<SplinepyBase>
  SplinepyExtractDim(const int& phys_dim) const;

  /// Derivative of composition
  virtual std::shared_ptr<SplinepyBase> SplinepyCompositionDerivative(
      const std::shared_ptr<SplinepyBase>& inner,
      const std::shared_ptr<SplinepyBase>& inner_derivative) const;

  /// @brief Gets proximity
  constexpr Proximity_& GetProximity() { return *proximity_; }
  /// @brief Gets proximity
  constexpr const Proximity_& GetProximity() const { return *proximity_; }

  /// Deep copy of current spline
  virtual std::shared_ptr<SplinepyBase> SplinepyDeepCopy() const {
    return std::make_shared<Bezier>(*this);
  };

protected:
  /// @brief Proximity
  std::shared_ptr<Proximity_> proximity_ = std::make_shared<Proximity_>(*this);
}; /* class Bezier */

} // namespace splinepy::splines

#include "splinepy/splines/bezier.inl"

#include "splinepy/explicit/bezier.hpp"
