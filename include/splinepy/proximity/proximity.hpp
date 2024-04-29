#pragma once

#include <napf.hpp>

#include "splinepy/splines/splinepy_base.hpp"
#include "splinepy/utils/arrays.hpp"
#include "splinepy/utils/grid_points.hpp"

namespace splinepy::proximity {

/// @brief compute convergence goal based on first norm and tolerances.
/// @tparam T
/// @param first_norm
/// @param rel_tol
/// @param abs_tol
/// @return
template<typename T>
T ConvergenceGoal(const T& first_norm, const T& rel_tol, const T& abs_tol) {
  return std::max(first_norm * rel_tol, abs_tol);
}

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
    int para_dim;
    int dim;

    struct {
      double J; // .5 distance^2
      double distance;
      double squared_distance;
      double convergence_norm;
      int stop_iteration;
    } status;

    struct {
      int max_iter;
      double tolerance;
      // this should be set by calling ComputeConvergenceGoal()
      double convergence_goal;
    } options;

    struct {
      double lambda;
      double prev_J;
      bool modified_update;
      double lower_bound;
      double upper_bound;
      RealArray_ candidate;

      /// @brief Computes || F(x) + J^T * dX ||^2
      /// @param spline_gradient
      /// @param delta_guess
      /// @param metric
      void Metric(const RealArray2D_& spline_gradient_,
                  const RealArray_& delta_guess_,
                  RealArray_& metric) const {
        for (int i{}; i < metric.size(); ++i) {
          metric[i] = delta_guess_.InnerProduct(&spline_gradient_(i, 0));
        }
      }
    } LM; // levenbergMarquart specific data

    struct {
      int m;      // n constraints
      int meq;    // n equality constraints
      int la;     // max(m, 1)
      int n;      // n variables
      double* x;  // current guess
      double* xl; // lower bounds
      double* xu; // upper bounds
      double f;   // value of objective
      double* c;  // m vector of c constraints
      double* g;  // partials of objective
      double* a;  // m x n constraint normals
      double acc; // abs(acc) is accuracy
      int iter;   // max iter
      int mode;   // mode -> see slsqp.f. 0 for entry
      double* w;  // work space -> see slsqp.f. need allocating
      int l_w;
      int* jw;
      int l_jw;
      double alpha;
      double f0;
      double gs;
      double h1;
      double h2;
      double h3;
      double h4;
      double t;
      double t0;
      double tol;
      int iexact;
      int incons;
      int ireset;
      int itermx;
      int line;
      int n1;
      int n2;
      int n3;

      RealArray_ w_arr;
      IndexArray_ jw_arr; // int array
      RealArray_ der_arr; //

      void SetWorkSpaces() {
        const int n1_ = n + 1;
        const int mineq = m - meq + n1 + n1;
        const int w_size = (3 * n1 + m) * (n1 + 1)
                           + (n1 - meq + 1) * (mineq + 2) + 2 * mineq
                           + (n1 + mineq) * (n1 - meq) + 2 * meq + n1
                           + (((n + 1) * n) / 2) + 2 * m + 3 * n + 3 * n1 + 1;

        // set w and lw
        w_arr.Reallocate(w_size);
        w = w_arr.data();
        l_w = w_size;

        // set jw and l_jw
        jw_arr.Reallocate(mineq);
        jw = jw_arr.data();
        l_jw = mineq;
      }

      void AllocateAndSetG() {
        der_arr.Reallocate(n + 1);
        der_arr[n - 1] = 0.0;
        g = der_arr.data();
      }
      void CopyDerivative(const RealArray_& der) {
        std::copy_n(der.data(), n, g);
      }
    } SLSQP; // slsqp specific data

    // view array
    ConstRealArray_ phys_query;
    // either view or allocated
    RealArray_ current_guess;
    RealArray_ current_phys;
    RealArray_ difference;
    RealArray2D_ spline_gradient;
    RealArray3D_ spline_hessian;

    // local allocated
    RealArray_ initial_guess;
    RealArray2D_ d2jdu2;
    RealArray_ djdu;
    RealArray_ delta_guess;
    RealArray2D_ spline_gradient_AAt;
    RealArray2D_ search_bounds;
    RealArray_ lower_bound; // view to search bounds
    RealArray_ upper_bound; // view to search bounds
    IndexArray_ clipped;

    /// lhs matrix based on 2D array. Created by search methods
    std::unique_ptr<LhsMatrix_> lhs_matrix;

    SearchData(const int para_dim,
               const int dim,
               const double tolerance,
               const int max_iter,
               const double* query,
               double* final_guess,
               double* nearest,
               double* nearest_minus_query,
               double* first_derivatives,
               double* second_derivatives) {
      Setup(para_dim,
            dim,
            tolerance,
            max_iter,
            query,
            final_guess,
            nearest,
            nearest_minus_query,
            first_derivatives,
            second_derivatives);
    }

    void Setup(const int para_dim_,
               const int dim_,
               const double tolerance,
               const int max_iter,
               const double* query,
               double* final_guess,
               double* nearest,
               double* nearest_minus_query,
               double* first_derivatives,
               double* second_derivatives) {

      para_dim = para_dim_;
      dim = dim_;

      options.tolerance = tolerance;
      options.max_iter = (max_iter < 0) ? para_dim * 20 : max_iter;
      options.convergence_goal = 0.0;

      // initialize values;
      status.J = -1.0;
      status.distance = -1.0;
      status.squared_distance = -1.0;
      status.convergence_norm = -1.0;

      phys_query.SetData(query, dim);

      current_guess.SetOrReallocateData(final_guess, para_dim);
      current_phys.SetOrReallocateData(nearest, dim);
      difference.SetOrReallocateData(nearest_minus_query, dim);
      spline_gradient.SetOrReallocateData(first_derivatives, para_dim, dim);
      spline_hessian.SetOrReallocateData(second_derivatives,
                                         para_dim,
                                         para_dim,
                                         dim);

      // allocate aux real arrays
      initial_guess.Reallocate(para_dim);
      d2jdu2.Reallocate(para_dim, para_dim);
      djdu.Reallocate(para_dim);
      delta_guess.Reallocate(para_dim);
      spline_gradient_AAt.Reallocate(para_dim, para_dim);
      search_bounds.Reallocate(2, para_dim);

      // get pointers to beginning of each bound
      lower_bound.SetData(search_bounds.begin(), para_dim);
      upper_bound.SetData(search_bounds.begin() + para_dim, para_dim);

      // allocate index arrays
      clipped.Reallocate(para_dim);
    }

    /// solvable lhs matrix - in mfem this is called Mult(rhs)
    LhsMatrix_& LhsMatrix() {
      if (!lhs_matrix) {
        splinepy::utils::PrintAndThrowError("lhs not set!");
      }
      return *lhs_matrix;
    }

    LhsMatrix_& LhsMatrix(RealArray2D_& arr) {
      lhs_matrix = std::make_unique<LhsMatrix_>(arr);
      return *lhs_matrix;
    }

    bool IsConverged() {
      return (status.convergence_norm < options.convergence_goal
              || status.distance < options.tolerance);
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

  /// @brief make initial guess - returns nearest neighbor from k-d tree
  /// @param goal
  /// @param guess
  void MakeInitialGuess(const ConstRealArray_& goal, RealArray_& guess) const;

  /// @brief make initial guess - returns nearest neighbor from k-d tree.
  /// Initial guess is saved to both initial_guess and current_guess
  /// @param goal
  /// @param guess
  void MakeInitialGuess(SearchData& aux) const;

  /// @brief updates current_guess with delta_guess and clips result
  /// @param aux
  void UpdateAndClip(SearchData& aux) const;

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

  /// @brief Fills dJ/du
  /// @param guess
  /// @param difference
  /// @param spline_gradient
  /// @param djdu
  void dJdu(const RealArray_& guess,
            const RealArray_& difference,
            RealArray2D_& spline_gradient,
            RealArray_& djdu) const;

  /// @brief Fills hessian of J
  /// @param guess
  /// @param difference
  /// @param spline_gradient_AAt
  /// @param spline_hessian
  /// @param lhs
  void d2Jdu2(const RealArray_& guess,
              const RealArray_& difference,
              const RealArray2D_& spline_gradient_AAt,
              RealArray3D_& spline_hessian,
              RealArray2D_& d2jdu2) const;

  /// @brief computes J and its derivatives based on the value of depth. depth
  /// == 1 -> jacobian, dept > 1 -> jacobian & hessian. Negative value is
  /// considered as highest depth we can compute (currently hessian)
  /// @param aux
  /// @param depth
  void ComputeCostAndDerivatives(SearchData& aux, int depth) const;

  /// @brief Computes status in SearchData based on current values in
  /// SearchData. If recompute_cose==true, it will compute the cost and its
  /// derivative again and store it.
  /// @param aux
  /// @param recompute_costs
  void ComputeStatus(SearchData& aux, const bool recompute_costs) const;

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

  /// @brief resets search bounds to spline's parametric bounds. if tight==true,
  /// current_guess +- grid_points_'s sampling step size.
  /// @param aux
  void FindSearchBound(SearchData& aux, const bool tight = false) const;

  /// @brief Prepares lhs and rhs of newton iteration
  void PrepareIterationNewton(SearchData& aux) const;

  /// @brief Performs newton iterations using PrepareIterationNewton
  /// @param aux
  void Newton(SearchData& aux) const;

  /// @brief Prepares lhs and rhs of LM iteration
  /// @param aux
  void PrepareIterationLevenbergMarquart(SearchData& aux) const;

  /// @brief Levenberg-Marquart method. Intended to be used as a fallback method
  /// for newton. Uses PrepareIterationLevenbergMarquart()
  /// @param aux
  void LevenbergMarquart(SearchData& aux) const;

  /// @brief Prepares lhs and rhs of slsqp iteration and other parameters
  void PrepareIterationSlsqp(SearchData& aux) const;

  /// @brief Performs slsqp iterations using PrepareIterationSlsqp
  /// @param aux
  void Slsqp(SearchData& aux) const;

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

/// @brief helper function to compute hessian at the end
/// @param spline
/// @param at
/// @param hess
void FillHessian(const splinepy::splines::SplinepyBase& spline,
                 const Proximity::RealArray_& at,
                 Proximity::RealArray3D_& hess);

} // namespace splinepy::proximity
