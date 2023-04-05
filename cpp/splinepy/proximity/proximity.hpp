#pragma once

#include <napf.hpp>

#include <splinepy/splines/helpers/properties.hpp>
#include <splinepy/utils/arrays.hpp>
#include <splinepy/utils/grid_points.hpp>
#include <splinepy/utils/nthreads.hpp>
#include <splinepy/utils/print.hpp>

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
template<typename SplineType>
class Proximity {
public:
  /// Options for initial guess
  enum class InitialGuess : int { MidPoint = 0, KdTree = 1 };

  // Frequently used array alias
  using DArrayD_ = std::array<double, SplineType::kDim>;
  using PArrayD_ = std::array<double, SplineType::kParaDim>;
  using DArrayI_ = std::array<int, SplineType::kDim>;
  using PArrayI_ = std::array<int, SplineType::kParaDim>;
  using PxPMatrixD_ = std::array<std::array<double, SplineType::kParaDim>,
                                 SplineType::kParaDim>;
  using PxDMatrixD_ =
      std::array<std::array<double, SplineType::kDim>, SplineType::kParaDim>;
  // kdtree related alias
  // C-Array instead of vector to avoid default init
  using Coordinates_ = typename SplineType::Coordinate_[];
  using Cloud_ = typename napf::
      CoordinatesCloud<std::unique_ptr<Coordinates_>, int, SplineType::kDim>;
  // metric is L2 and returned distance will be squared.
  using Tree_ =
      typename std::conditional<(SplineType::kDim < 4),
                                napf::CoordinatesTree<double, /* DataT */
                                                      double, /* DistT */
                                                      int,    /* IndexT */
                                                      SplineType::kDim, /* dim
                                                                         */
                                                      2,       /* metric (L2) */
                                                      Cloud_>, /* Cloud T */
                                napf::CoordinatesHighDimTree<double,
                                                             double,
                                                             int,
                                                             SplineType::kDim,
                                                             2,
                                                             Cloud_>>::type;
  using GridPoints_ =
      splinepy::utils::GridPoints<double, int, SplineType::kParaDim>;

  /// Constructor. As a spline helper class, always need a spline.
  Proximity(SplineType const& spline) : spline_(spline){};

  /*!
   * Computes difference between physical query and current parametric guess.
   *
   * @param[in] query
   * @param[in] guess
   * @param[out] difference query - guess
   */
  void GuessMinusQuery(const typename SplineType::ParametricCoordinate_& guess,
                       const double* query,
                       DArrayD_& difference) const {

    splinepy::utils::FirstMinusSecondEqualsThird(spline_(guess),
                                                 query,
                                                 difference);
  }

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
  void PlantNewKdTree(const PArrayI_& resolutions, const int n_thread = 1) {
    // skip early, if requested resolutions are the as existing kdtree
    if (kdtree_planted_ && sampled_resolutions_ == resolutions) {
      return;
    }

    // create fresh grid_points_ and coordinates_
    const auto parametric_bounds =
        splinepy::splines::helpers::GetParametricBounds(spline_);
    grid_points_ = GridPoints_(parametric_bounds, resolutions);
    coordinates_ = std::make_unique<Coordinates_>(grid_points_.Size());

    // lambda function to allow n-thread execution
    auto sample_coordinates = [&](int begin, int end) {
      typename SplineType::ParametricCoordinate_ parametric_coordinate;
      for (int i{begin}; i < end; ++i) {
        grid_points_.IndexToGridPoint(i, parametric_coordinate);
        coordinates_[i] = spline_(parametric_coordinate);
      }
    };

    // n-thread execution
    splinepy::utils::NThreadExecution(sample_coordinates,
                                      grid_points_.Size(),
                                      n_thread);

    // plant a new tree
    coordinates_cloud_ =
        std::make_unique<Cloud_>(coordinates_, grid_points_.Size());
    kdtree_ = std::make_unique<Tree_>(SplineType::kDim, *coordinates_cloud_);
    kdtree_planted_ = true;
  }

  /// Make initial guess of choice
  typename SplineType::ParametricCoordinate_
  MakeInitialGuess(const int& initial_guess, const double* goal) const {

    if (initial_guess == static_cast<int>(InitialGuess::MidPoint)) {
      const auto parametric_bounds =
          splinepy::splines::helpers::GetParametricBounds(spline_);

      typename SplineType::ParametricCoordinate_ mean;
      splinepy::utils::Mean(parametric_bounds[0], parametric_bounds[1], mean);

      return mean;

    } else if (initial_guess == static_cast<int>(InitialGuess::KdTree)) {
      if (!kdtree_planted_) {
        // hate to be aggresive, but here it is.
        splinepy::utils::PrintAndThrowError(
            "to use InitialGuess::Kdtree option,"
            "please first plant a kdtree.",
            "For example:\n",
            "  SplineType spline{ ... /* spline init */ };\n",
            "  std::array<int, SplineType::kParaDim>",
            "resolutions{ ... /* kdtree sample resolutions*/ };\n",
            "  const int nthreads = ... /* number of threads */;\n",
            "  spline.GetProximity().PlantNewKdTree(resolutions, nthreads);\n",
            "\n  /* For SplinepyBase */\n"
            "  SplinepyBase spline{... /* splinepybase init */};\n",
            "  std::vector<int> resolutions(spline.SplinepyParaDim());\n",
            "  ... /* fill resolutions */ ...\n",
            "  const int nthreads = ... /* number of threads */;\n",
            "  spline.SplinepyPlantNewKdtreeForProximity(resolutions.data(),",
            "nthreads);\n");
      }

      // good to go. ask the tree
      int id;
      double distance;
      kdtree_->knnSearch(goal, 1 /* closest neighbor */, &id, &distance);

      using ReturnType = typename SplineType::ParametricCoordinate_;

      return grid_points_.template IndexToGridPoint<ReturnType>(id);
    } else {
      // well, shouldn't reach here, but to fight a warning, here it is
      splinepy::utils::PrintAndThrowError("Invalid option for initial guess!");
      return typename SplineType::ParametricCoordinate_{};
    }
  }

  /*!
   * Builds RHS and fills spline_gradient, which is also required in LHS.
   *
   * @params[in] guess current parametric coordinate guess
   * @params[in] difference result of `GuessMinusQuery()`
   * @params[out] spline_gradient
   * @params[out] rhs
   */
  void FillSplineGradientAndRhs(
      const typename SplineType::ParametricCoordinate_& guess,
      const DArrayD_& difference,
      PxDMatrixD_& spline_gradient,
      PArrayD_& rhs) const {
    // btw, RHS is what's internally called as "df_dxi"
    double df_dxi_i; // temporary variable to hold sum

    typename SplineType::Derivative_ derivative_query;
    using DerivativeValueType = typename SplineType::Derivative_::value_type;

    for (int i{}; i < SplineType::kParaDim; ++i) {
      // spline derivative query formulation
      derivative_query.fill(DerivativeValueType{0});
      ++derivative_query[i];
      // derivative evaluation
      auto const derivative = spline_(guess, derivative_query);

      df_dxi_i = 0.;
      for (int j{}; j < SplineType::kDim; ++j) {
        // cast, since d may be a NamedType.
        double d_value;
        if constexpr (std::is_scalar<decltype(derivative)>::value) {
          d_value = derivative;
        } else {
          d_value = static_cast<double>(derivative[j]);
        }
        df_dxi_i += difference[j] * d_value;
        spline_gradient[i][j] = d_value;
      }
      rhs[i] = -2. * df_dxi_i; // apply minus here already!
    }
  }

  /*!
   * Builds LHS
   *
   * @params[in] guess
   * @params[in] difference
   * @params[in] spline_gradient
   * @params[out] lhs
   */
  void FillLhs(const typename SplineType::ParametricCoordinate_& guess,
               const DArrayD_& difference,
               const PxDMatrixD_& spline_gradient,
               PxPMatrixD_& lhs) const {

    const PxPMatrixD_ spline_gradientAAt =
        splinepy::utils::AAt(spline_gradient);

    typename SplineType::Derivative_ derivative_query;
    using DerivativeValueType = typename SplineType::Derivative_::value_type;
    for (int i{}; i < SplineType::kParaDim; ++i) {
      for (int j{i}; j < SplineType::kParaDim; ++j) {
        derivative_query.fill(DerivativeValueType{0});
        ++derivative_query[i];
        ++derivative_query[j];

        auto const derivative = spline_(guess, derivative_query);

        lhs[i][j] = 2
                    * (splinepy::utils::Dot(difference, derivative)
                       + spline_gradientAAt[i][j]);

        if (i != j)
          lhs[j][i] = lhs[i][j]; // fill symmetric part
      }
    }
  }

  /*!
   * Given physical coordinate, finds closest parametric coordinate.
   *
   * @params query
   * @params initial_guess
   * @params tolerance
   */
  typename SplineType::ParametricCoordinate_
  FindNearestParametricCoordinate(const double* query,
                                  InitialGuess initial_guess,
                                  double tolerance = 1e-12,
                                  bool aggressive_bounds = false) const {

    PxPMatrixD_ lhs;
    PArrayD_ rhs;
    PArrayD_ delta_guess;
    DArrayD_ difference;
    PxDMatrixD_ spline_gradient;
    PArrayI_ clipped{};          /* clip status after most recent update */
    PArrayI_ previous_clipped{};
    PArrayI_ solver_skip_mask{}; /* tell solver to skip certain entry */
    bool solver_skip_mask_activated = false;

    // search_bounds is parametric bounds.
    auto search_bounds =
        splinepy::splines::helpers::GetParametricBounds(spline_);

    typename SplineType::ParametricCoordinate_ current_guess =
        MakeInitialGuess(static_cast<int>(initial_guess), query);

    // Be optimistic and check if initial guess was awesome
    // If it is awesome, return and have feierabend
    GuessMinusQuery(current_guess, query, difference);
    if (splinepy::utils::NormL2(difference) < tolerance) {
      return current_guess;
    }

    // Let's try aggresive search bounds
    if (initial_guess == InitialGuess::KdTree && aggressive_bounds) {
      // you need to be sure that you have sampled your spline fine enough
      for (std::size_t i{}; i < SplineType::kParaDim; ++i) {
        // adjust lower (0) and upper (1) bounds aggressively
        // but of course, not so aggresive that it is out of bound.
        search_bounds[0][i] =
            std::max(search_bounds[0][i],
                     current_guess[i] - grid_points_.step_size_[i]);
        search_bounds[1][i] =
            std::min(search_bounds[1][i],
                     current_guess[i] + grid_points_.step_size_[i]);
      }
    }

    const int max_iteration = SplineType::kParaDim * 20;
    double previous_norm{}, current_norm;

    // newton iteration
    for (int i{}; i < max_iteration; ++i) {
      // build systems to solve
      FillSplineGradientAndRhs(current_guess, difference, spline_gradient, rhs);
      FillLhs(current_guess, difference, spline_gradient, lhs);

      // solve and update
      // 1. set solver skip mask if clipping happend twice at the same place.
      if (previous_clipped == clipped && !solver_skip_mask_activated) {
        solver_skip_mask = clipped;
        solver_skip_mask_activated = true;
        // if skip mask is on for all entries, return now.
        if (splinepy::utils::NonZeros(solver_skip_mask)
            == SplineType::kParaDim) {
          // current_guess should be clipped at this point.
          return current_guess;
        }
      }

      // GaussWithPivot may swap and modify enties of all the input
      // -> can't use lhs and rhs afterwards, and we don't need them.
      // -> solver_skip_mask and delta_guess is reordered to rewind swaps
      splinepy::utils::GaussWithPivot(lhs, rhs, solver_skip_mask, delta_guess);
      // Update
      splinepy::utils::AddSecondToFirst(current_guess, delta_guess);
      // Clip
      previous_clipped = clipped; // swap?
      splinepy::utils::Clip(search_bounds, current_guess, clipped);
      // Converged?
      current_norm = splinepy::utils::NormL2(rhs);
      GuessMinusQuery(current_guess, query, difference);
      if (std::abs(previous_norm - current_norm) < tolerance
          || splinepy::utils::NormL2(difference) < tolerance)
        break;

      // prepare next round
      previous_norm = current_norm;
    }
    return current_guess;
  }

  void FirstOrderFallBack() {}

  /*!
   * Given physical coordinate, finds closest parametric coordinate.
   * Always takes initial guess based on kdtree.
   *
   * @params[in] query
   * @params[in] tolerance
   * @params[in] max_iterations
   * @params[in] aggresive_bounds
   * @params[out] final_guess (dim)
   * @params[out] nearest (dim)
   * @params[out] nearest_minus_query (dim)
   * @params[out] distance
   * @params[out] convergence_norm
   * @params[out] first_derivatives (para_dim x dim)
   * @params[out] second_derivatives (para_dim x para_dim x dim)
   */
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
                    double* second_derivatives /* spline hessian */) const {

    PxPMatrixD_ lhs;
    PArrayD_ rhs;
    PArrayD_ delta_guess;
    DArrayD_ difference;
    PxDMatrixD_ spline_gradient;
    PArrayI_ clipped{};          /* clip status after most recent update */
    PArrayI_ previous_clipped{};
    PArrayI_ solver_skip_mask{}; /* tell solver to skip certain entry */
    typename SplineType::Coordinate_ current_phys{};
    double current_distance{};
    double previous_norm{}, current_norm{};

    // search_bounds is parametric bounds here
    auto search_bounds =
        splinepy::splines::helpers::GetParametricBounds(spline_);

    // for verbose, we don't return right away even this is the best guess
    // already, so that we can fill out all the other infos.
    typename SplineType::ParametricCoordinate_ current_guess =
        MakeInitialGuess(static_cast<int>(InitialGuess::KdTree), query);

    // Let's try aggresive search bounds
    if (aggressive_bounds) {
      // you need to be sure that you have sampled your spline fine enough
      for (std::size_t i{}; i < SplineType::kParaDim; ++i) {
        // adjust lower (0) and upper (1) bounds aggressively
        // but of course, not so aggresive that it is out of bound.
        search_bounds[0][i] =
            std::max(search_bounds[0][i],
                     current_guess[i] - grid_points_.step_size_[i]);
        search_bounds[1][i] =
            std::min(search_bounds[1][i],
                     current_guess[i] + grid_points_.step_size_[i]);
      }
    }
    const int max_iter =
        max_iterations < 0 ? SplineType::kParaDim * 20 : max_iterations;

    // build systems to solve
    GuessMinusQuery(current_guess, query, difference);
    FillSplineGradientAndRhs(current_guess, difference, spline_gradient, rhs);

    // 0 iteration returns initial guess.
    // compute rest of verbose info here.
    if (max_iterations == 0) {
      current_phys = spline_(current_guess);
      splinepy::utils::FirstMinusSecondEqualsThird(current_phys,
                                                   query,
                                                   difference);
      current_distance = splinepy::utils::NormL2(difference);
      current_norm = std::abs(splinepy::utils::NormL2(rhs));
    }

    // newton iterations
    for (int i{}; i < max_iter; ++i) {
      // lhs
      FillLhs(current_guess, difference, spline_gradient, lhs);

      // GaussWithPivot may swap and modify enties of all the input
      // -> can't use lhs and rhs afterwards, and we don't need them.
      // -> solver_skip_mask and delta_guess is reordered to rewind swaps
      splinepy::utils::GaussWithPivot(lhs, rhs, solver_skip_mask, delta_guess);
      // Update
      splinepy::utils::AddSecondToFirst(current_guess, delta_guess);
      // Clip
      splinepy::utils::Clip(search_bounds, current_guess, clipped);
      // check distance
      current_phys = spline_(current_guess);
      splinepy::utils::FirstMinusSecondEqualsThird(current_phys,
                                                   query,
                                                   difference);
      current_distance = splinepy::utils::NormL2(difference);
      // assemble rhs, check norm
      FillSplineGradientAndRhs(current_guess, difference, spline_gradient, rhs);
      current_norm = splinepy::utils::NormL2(rhs);

      // convergence check
      if (std::abs(previous_norm - current_norm) < tolerance
          || current_distance < tolerance) {
        break;
      }
      // set solver skip mask if clipping happened twice at the same place.
      if (previous_clipped == clipped) {
        solver_skip_mask = clipped;
        // if skip mask is on for all entries, return now.
        if (splinepy::utils::NonZeros(solver_skip_mask)
            == SplineType::kParaDim) {
          // current_guess should be clipped at this point.
          break;
        }
      }

      // we are here because it didn't converge. prepare next round
      previous_norm = current_norm;
      std::swap(previous_clipped, clipped);
    }
    // write return values - 7 args
    distance = current_distance;     /* 1 */
    convergence_norm = current_norm; /* 2 */
    typename SplineType::Derivative_ derivative_query;
    double der; /* to accomodate different kind of derivative types */
    using DerivativeValueType = typename SplineType::Derivative_::value_type;
    for (int i{}; i < SplineType::kParaDim; ++i) {
      final_guess[i] = static_cast<double>(current_guess[i]); /* 3 */
      for (int j{i}; j < SplineType::kParaDim; ++j) {
        derivative_query.fill(DerivativeValueType{0});
        ++derivative_query[i];
        ++derivative_query[j];
        const auto derivative = spline_(current_guess, derivative_query);
        for (int k{}; k < SplineType::kDim; ++k) {
          if constexpr (std::is_scalar<decltype(derivative)>::value) {
            der = derivative;
          } else {
            der = static_cast<double>(derivative[k]);
          }
          // spline hessian
          second_derivatives[(i * SplineType::kParaDim * SplineType::kDim)
                             + (j * SplineType::kDim) + k] = der; /* 4 */
          // symmetric part
          if (i != j) {
            second_derivatives[(j * SplineType::kParaDim * SplineType::kDim)
                               + (i * SplineType::kDim) + k] = der;
          }

          // ones that don't need extra para_dim loop
          if (i == 0 /* j starts with 0 */) {
            first_derivatives[j * SplineType::kDim + k] =
                spline_gradient[j][k]; /* 5 */
            // ones that don't need extra extra para_dim loop
            // => pure dim loop
            if (j == 0) {
              double cur_phys;
              if constexpr (std::is_scalar_v<decltype(current_phys)>) {
                nearest[k] = current_phys;
              } else {
                nearest[k] = static_cast<double>(current_phys[k]); /* 6 */
              }
              nearest_minus_query[k] = difference[k];              /* 7 */
            }
          }
        }
      }
    }
  }

protected:
  SplineType const& spline_;
  // kdtree related variables
  GridPoints_ grid_points_;
  PArrayI_ sampled_resolutions_;
  bool kdtree_planted_ = false;
  std::unique_ptr<Coordinates_> coordinates_;
  std::unique_ptr<Cloud_> coordinates_cloud_;
  std::unique_ptr<Tree_> kdtree_;
};

} // namespace splinepy::proximity
