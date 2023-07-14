#pragma once

#include <array>
#include <type_traits>

// bezman
#include <bezman/src/bezier_group.hpp>
#include <bezman/src/bezier_spline.hpp>
#include <bezman/src/point.hpp>
#include <bezman/src/vertices.hpp>

#include <splinepy/explicit/bezman/bezier_extern.hpp>

#include <splinepy/proximity/proximity.hpp>
#include <splinepy/splines/helpers/basis_functions.hpp>
#include <splinepy/splines/helpers/extract.hpp>
#include <splinepy/splines/helpers/properties.hpp>
#include <splinepy/splines/helpers/scalar_type_wrapper.hpp>
#include <splinepy/splines/splinepy_base.hpp>
#include <splinepy/utils/print.hpp>

namespace splinepy::splines {

template<std::size_t para_dim>
using BezierSplineType =
    bezman::BezierSpline<para_dim, bezman::VertexView<double>, double>;

// splinepy's Rational Bezier
/// @brief Rational Bezier class
/// @tparam para_dim Dimension of parametric space
/// @tparam dim Dimension of physical space
//template<std::size_t para_dim>
//class RationalBezier;

/// @brief Bezier class
/// @tparam para_dim Dimension of parametric space
/// @tparam dim Dimension of physical space
template<std::size_t para_dim>
class Bezier : public SplinepyBase,
               public BezierSplineType<para_dim> {
public:
  /// @brief Parametric space dimension
  static constexpr int kParaDim = static_cast<int>(para_dim);
  /// @brief True iff it is rational Bezier
  static constexpr bool kIsRational = false;
  /// @brief True iff spline has knot vectors
  static constexpr bool kHasKnotVectors = false;

  using SplinepyBase_ = SplinepyBase;
  template<std::size_t s_para_dim>
  using SelfTemplate_ = Bezier<s_para_dim>;
  using SelfBoundary_ = Bezier<para_dim - 1>;
  using Base_ = BezierSplineType<para_dim>;
  template<std::size_t b_para_dim>
  using BaseTemplate_ = BezierSplineType<b_para_dim>;
  // alias to enable helper functions.
  using ParametricCoordinate_ = typename bezman::Point<para_dim, double>;
  using Degrees_ = typename std::array<std::size_t, para_dim>;
  // using Coordinate_ = typename Base_::PhysicalPointType_;
  //  keeping the notation similar to BSplineLib
  using VectorSpace_ = typename Base_::ControlPointsType_;
  using Coordinate_ = typename VectorSpace_::ContainerType_;
  using Coordinates_ = typename VectorSpace_::ContainerType_; // flat vector
  using Derivative_ = typename std::array<std::size_t, para_dim>;
  using Dimension_ = std::size_t;
  using ProximityBase_ = proximity::ProximityBase;
  template<int dim>
  using Proximity_ = proximity::Proximity<Bezier<para_dim>, dim>;

  /// @brief Creates Base for Bezier
  /// @param degrees
  /// @param control_points
  static Base_
  CreateBase(const int* degrees, const double* control_points, const int dim);

  /// @brief Constructor based on raw pointer
  /// @param degrees Degrees of the Bezier
  /// @param control_points Control points of the Bezier
  Bezier(const int* degrees, const double* control_points, const int dim)
      : Base_(CreateBase(degrees, control_points, dim)) {}

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
  virtual int SplinepyDim() const {
    return Base_::control_points.GetDimension();
  }

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
  virtual int SplinepyNumberOfControlPoints() const {
    return static_cast<int>(Base_::control_points.size());
  }

  /// @copydoc splinepy::splines::SplinepyBase::SplinepyNumberOfSupports
  virtual int SplinepyNumberOfSupports() const {
    return splinepy::splines::helpers::GetNumberOfSupports(*this);
  }

  /// @copydoc splinepy::splines::SplinepyBase::SplinepyCurrentProperties
  virtual void SplinepyCurrentProperties(
      int* degrees,
      std::vector<std::vector<double>>* knot_vectors /* untouched */,
      double* control_points,
      double* weights /* untouched */) const;

  virtual void SplinepyParametricBounds(double* para_bounds) const;

  virtual void SplinepyControlMeshResolutions(int* control_mesh_res) const;

  /// @brief Calculate Greville abscissae for Bezier
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

  constexpr void CreateProximityIfNull();

  /// @brief Gets proximity
  constexpr ProximityBase_& GetProximity();

  /// @brief Gets proximity
  constexpr const ProximityBase_& GetProximity() const;

protected:
  /// @brief Proximity
  std::shared_ptr<ProximityBase_> proximity_ = nullptr;

}; /* class Bezier */

//#include <splinepy/splines/bezier.inl>

} // namespace splinepy::splines

#include <splinepy/explicit/splinepy/bezier_extern.hpp>
