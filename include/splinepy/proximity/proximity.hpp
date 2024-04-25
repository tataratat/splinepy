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
 *
 * J = 1 / 2 * (p - S(u))^2
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

  using LhsMatrix_ = splinepy::utils::Matrix<double, int>;

  /// struct to hold temporary data during the search.
  /// This should allow easier transition to fallback methods.
  struct SearchData {
    // view array
    ConstRealArray_ phys_query_;
    // either view or allocated
    RealArray_ current_guess_;
    RealArray_ current_phys_;
    RealArray_ difference_;
    RealArray2D_ spline_gradient_;
    RealArray3D_ spline_hessian_;

    // local allocated
    RealArray_ initial_guess_;
    RealArray2D_ lhs_;
    RealArray_ rhs_;
    RealArray_ delta_guess_;
    RealArray2D_ spline_gradient_AAt_;
    RealArray2D_ search_bounds_;
    RealArray_ lower_bound_; // view to search bounds
    RealArray_ upper_bound_; // view to search bounds
    IndexArray_ clipped_;

    std::unique_ptr<LhsMatrix_> lhs_matrix_;

    SearchData(const int para_dim,
               const int dim,
               const double* query,
               double* final_guess,
               double* nearest,
               double* nearest_minus_query,
               double* first_derivatives,
               double* second_derivatives) {
      Setup(para_dim,
            dim,
            query,
            final_guess,
            nearest,
            nearest_minus_query,
            first_derivatives,
            second_derivatives);
    }

    void Setup(const int para_dim,
               const int dim,
               const double* query,
               double* final_guess,
               double* nearest,
               double* nearest_minus_query,
               double* first_derivatives,
               double* second_derivatives) {

      phys_query_.SetData(query, dim);

      current_guess_.SetOrReallocateData(final_guess, para_dim);
      current_phys_.SetOrReallocateData(nearest, dim);
      difference_.SetOrReallocateData(nearest_minus_query, dim);
      spline_gradient_.SetOrReallocateData(first_derivatives, para_dim, dim);
      spline_hessian_.SetOrReallocateData(second_derivatives,
                                          para_dim,
                                          para_dim,
                                          dim);

      // allocate aux real arrays
      initial_guess_.Reallocate(para_dim);
      lhs_.Reallocate(para_dim, para_dim);
      rhs_.Reallocate(para_dim);
      delta_guess_.Reallocate(para_dim);
      spline_gradient_AAt_.Reallocate(para_dim, para_dim);
      search_bounds_.Reallocate(2, para_dim);

      // get pointers to beginning of each bound
      lower_bound_.SetData(search_bounds.begin(), para_dim);
      upper_bound_.SetData(search_bounds.begin() + para_dim, para_dim);

      // allocate index arrays
      clipped_.Reallocate(para_dim);
    }

    /// solvable lhs matrix - in mfem this is called Mult(rhs)
    LhsMatrix_& LhsMatrix() {
      if (!lhs_matrix_) {
        lhs_matrix_ = std::make_unique<LhsMatrix_>(lhs_);
      }
      return *lhs_matrix_;
    }
  };

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

  /// @brief returns J = .5 * distance^2. Where distance is norm of difference.
  /// Difference is defined as || query - spline(guess_phys) ||
  /// @param guess
  /// @param query
  /// @param guess_phys
  /// @param difference
  /// @return
  double J(const RealArray_& guess,
           const ConstRealArray_& query,
           RealArray_& guess_phys,
           RealArray_& difference) const;

  /// @brief returns J = .5 * distance^2. Where distance is norm of difference.
  /// Difference is defined as || query - spline(guess_phys) ||
  /// @param aux
  /// @return
  double J(SearchData& aux) const;

  /// @brief Fills dJ/du = difference * dS/du
  /// @param guess
  /// @param difference
  /// @param spline_gradient
  /// @param djdu
  void dJdu(const RealArray_& guess,
            const RealArray_& difference,
            RealArray2D_& spline_gradient,
            RealArray2D_& djdu) const;

  /// @brief Fills dJ/du = difference * dS/du
  void dJdu(SearchData& aux) const;

  /// @brief Fills hessian of J d2J/dujdjk = dS/duk dS/duj +
  /// difference*dS/dujduk
  /// @param guess
  /// @param difference
  /// @param spline_gradient_AAt
  /// @param spline_hessian
  /// @param lhs
  void d2Jdu2(const RealArray_& guess,
              const RealArray_& difference,
              const RealArray2D_& spline_gradient_AAt,
              RealArray3D_& spline_hessian,
              RealArray2D_& lhs) const;

  /// @brief Fills hessian of J d2J/dujdjk = dS/duk dS/duj +
  /// difference*dS/dujduk
  /// @param aux
  void d2Jdu2(SearchData& aux) const;

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

  /// Fill LHS of Newton method for critical value search
  /// Note that only upper triangle is filled for spline hessian
  void FillLhs(const RealArray_& guess,
               const RealArray_& difference,
               const RealArray2D_& spline_gradient_AAt,
               RealArray3D_& spline_hessian,
               RealArray2D_& lhs) const;

  // Fill LHS of damped Levenberg-Marquart method
  void FillLhsLM(const RealArray_& guess,
                 const RealArray_& difference,
                 const RealArray2D_& spline_gradient_AAt,
                 RealArray2D_& lhs,
                 const double& lambda,
                 const bool& modified_marquart) const;

  /// @brief Prepares lhs and rhs of LM iteration
  /// @param aux
  void PrepareIterationLevenbergMarquart(SearchData& aux) const;

  /// @brief Levenberg-Marquart method. Intended to be used as a fallback method
  /// for newton. Uses PrepareIterationLevenbergMarquart()
  /// @param aux
  void LevenbergMarquart(SearchData& aux) const;

  /// @brief Prepares lhs and rhs of newton iteration
  void PrepareIterationNewton(SearchData& aux) const;

  /// @brief Performs newton iterations using PrepareIterationNewton
  /// @param aux
  void Newton(SearchData& aux) const;

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
