#pragma once

#include <napf.hpp>

#include "splinepy/splines/splinepy_base.hpp"
#include "splinepy/utils/arrays.hpp"
#include "splinepy/utils/grid_points.hpp"

namespace splinepy::proximity {

/*!
 * A helper class to perform proximity operations for splines.
 *
 * Given a physical coordinate, this tries to find parametric coordinate that
 * maps to the closest physical coordinate. Often referred as
 * "point inversion". This is done by searching for a root of first derivative
 * of squared distance between mapped coordinate and query coordinate.
 * For detailed information, please take a look at splinepy python
 * documentation.
 */
class Proximity {
public:
  /// @brief array cloud that wraps data of Array
  using Cloud_ = napf::ArrayCloud<double, int>;
  /// @brief metric is L2 and this returns squared distance
  using Tree_ = napf::ArrayTree<double, double, int, 2>;

  using RealArray_ = splinepy::utils::Array<double>;
  using RealArray2D_ = splinepy::utils::Array<double, 2>;
  using RealArray3D_ = splinepy::utils::Array<double, 3>;
  using IndexArray_ = splinepy::utils::Array<int>;

  // for non-writable views
  using ConstRealArray_ = splinepy::utils::Array<const double>;

  using SystemMatrix = splinepy::utils::Matrix<double, int>;

protected:
  // helpee spline
  const splinepy::splines::SplinepyBase& spline_;

  // kdtree related variables
  splinepy::utils::GridPoints grid_points_;
  RealArray_ sampled_spline_;
  std::unique_ptr<Cloud_> cloud_;
  std::unique_ptr<Tree_> kdtree_;

public:
  /// Constructor. As a spline helper class, always need a spline.
  Proximity(const splinepy::splines::SplinepyBase& spline) : spline_(spline){};

  /*!
   * Plants a kdtree with given resolution.
   *
   * This needs to be built before BEFORE you request a proximity query with
   * `InitialGuess::kdTree`. This will always plant a new tree: at runtime, if
   * a finer tree is desired, you can plant it again.
   *
   * @param resolutions parameter space sampling resolution
   * @param n_thread number of threads to be used for sampling
   */
  void PlantNewKdTree(const int* resolutions, const int n_thread = 1);

  /// @brief difference = spline(guess) - query. In current formulation, this is
  /// our objective function.
  /// @param guess
  /// @param query
  /// @param difference
  void GuessMinusQuery(const RealArray_& guess,
                       const ConstRealArray_& query,
                       RealArray_& difference) const;

  /// @brief difference = spline(guess) - query.
  /// This returns spline(guess)
  /// @param guess
  /// @param query
  /// @param guess_phys
  /// @param difference
  void GuessMinusQuery(const RealArray_& guess,
                       const ConstRealArray_& query,
                       RealArray_& guess_phys,
                       RealArray_& difference) const;

  void MakeInitialGuess(const ConstRealArray_& goal, RealArray_& guess) const;

  /*!
   * Builds RHS and fills spline_gradient, which is also required in LHS.
   * RHS is what's internally called as "df_dxi"
   *
   * @param[in] guess current parametric coordinate guess
   * @param[in] difference result of `GuessMinusQuery()`
   * @param[out] spline_gradient
   * @param[out] rhs
   */
  void FillSplineGradientAndRhs(const RealArray_& guess,
                                const RealArray_& difference,
                                RealArray2D_& spline_gradient,
                                RealArray_& rhs) const;

  void FillLhs(const RealArray_& guess,
               const RealArray_& difference,
               const RealArray2D_& spline_gradient_AAt,
               RealArray3D_& spline_hessian,
               RealArray2D_& lhs) const;

  /// @brief First order fall back
  void FirstOrderFallBack() {}

  /// @brief Given physical coordinate, finds closest parametric coordinate.
  /// Always takes initial guess based on kdtree.
  ///
  /// @param[in] query
  /// @param[in] tolerance
  /// @param[in] max_iterations
  /// @param[in] aggressive_bounds
  /// @param[out] final_guess (para_dim)
  /// @param[out] nearest (dim)
  /// @param[out] nearest_minus_query (dim)
  /// @param[out] distance
  /// @param[out] convergence_norm
  /// @param[out] first_derivatives (para_dim x dim)
  /// @param[out] second_derivatives (para_dim x para_dim x dim)
  void VerboseQuery(const double* query,
                    const double& tolerance,
                    const int& max_iterations,
                    const bool aggressive_bounds,
                    double* final_guess,
                    double* nearest /* spline(final_guess) */,
                    double* nearest_minus_query /* difference */,
                    double& distance,
                    double& convergence_norm,
                    double* first_derivatives /* spline jacobian */,
                    double* second_derivatives /* spline hessian */) const;
};

} // namespace splinepy::proximity
