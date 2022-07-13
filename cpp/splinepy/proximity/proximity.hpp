#pragma once

#include <utils/arrays.hpp>
#include <utils/print.hpp>
#include <utils/nthreads.hpp>
#include <utils/grid_points.hpp>

namespace splinepy::proximity {


/*!
 * A helper class to perform proximity operations for splines.
 */
template<typename SplineType>
class Proximity {
public:

  /// Options for initial guess
  enum class InitalGuess {MidPoint, KdTree};

  // Frequently used array alias
  using DArrayD_ = std::array<double, SplineType::kDim>;
  using PArrayD_ = std::array<double, SplineType::kParaDim>;
  using DArrayI_ = std::array<int, SplineType::kDim>;
  using PArrayI_ = std::array<int, SplineType::kParaDim>;
  using PxPMatrixD_ = std::array<
      std::array<double, SplineType::kParaDim>, SplineType::kParaDim
  >;
  using PxDMatrixD_ = std::array<
      std::array<double, SplineType::kDim>, SplineType::kParaDim
  >;
  // kdtree related alias
  using Coordinates_ = typename std::unique_ptr<SplineType::Coordinate_[]>;
  using Cloud_ =
      typename napf::CoordinatesCloud<Coordinates_, int, SplineType::kDim>;
  using Tree_ = typename std::conditional<
    (SplineType::kDim < 4),
    napf::CoordinatesTree<double,  /* DataT */
                          double,  /* DistT */
                          int,     /* IndexT */
                          dim,     /* dim */
                          1,       /* metric (L1) */
                          Cloud_>, /* CloudT */
    napf::CoordinatesHighDimTree<double, double, int, dim, 1, Cloud_>
  >::type;
  using GridPoints_ = GridPoints<double, int, SplineType::kParaDim>;

  /// Constructor. As a spline helper class, always need a spline.
  Proximity(SplineType const &spline) : spline_(spline) {};


  /*!
   * Computes difference between physical query and current parametric guess.
   *
   * @param[in] query
   * @param[in] guess
   * @param[out] difference query - guess
   */
  void QueryMinusGuess(const double* query,
                       const SplineType::ParametricCoordinate_& guess,
                       DArrayD_& difference) const {

    FirstMinusSecondEqualsThird(
        query,
        spline_.(guess),
        difference
    );

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
  void PlantNewKdTree(const PArrayI_& resolutions, const int n_thread=1) {
    // skip early, if requested resolutions are the as existing kdtree
    if (kdtree_platned_ && sampled_resolutions_ == resolutions) {
      return;
    }

    // create fresh grid_points_ and coordinates_
    const auto parametric_bounds =
        splinepy::splines::GetParametricBounds(splines_);
    grid_points_ = GridPoints_(parametric_bounds, resolutions);
    coordinates_ =
        Coordinates_(new SplineType::Coordinate_[grid_points_.Size()]);

    // lambda function to allow n-thread execution
    auto sample_coordinates = [&] (int begin, int end) {
      SplineType::ParametricCoordinate_ parametric_coordinate;
      for (int i{begin}; i < end; ++i) {
        grid_points_.IndexToParametricCoordinate(i, parametric_coordinate);
        coordinates_[i] = std::move(spline_(parametric_coordinate));
      }
    }

    // n-thread execution
    splinepy::utils::NThreadExecution(sample_coordinates,
                                      grid_points_.Size(),
                                      n_thread);

    // plant a new tree
    coordinates_cloud_ = std::make_unique<Cloud_>(coordinates_,
                                                  grid_points_.size());
    kdtree_ = std::make_unique<Tree_>(SplineType::kDim, *coordinates_cloud_);
    kdtree_planted_ = true;
  }

  /// Make initial guess of choice
  SplineType::ParametricCoordinate_ MakeInitialGuess(
      InitialGuess initial_guess,
      const double* goal) const {


    if (initial_guess == InitialGuess::MidPoint) {
      using ReturnValueType = SplineType::ParametricCoordinate_::value_type;

      // mid point is mean of parametric bounds. doesn't consider the goal.
      return splinepy::utils::Mean<ReturnValueType>(
          splinepy::splines::GetParametricBounds(splines_
      );

    } else (initial_guess == InitialGuess::KdTree) {
      if (!kdtree_planted_) {
        // hate to be aggresive, but here it is.
        splinepy::utils::PrintAndThrowError(
            "to use InitialGuess::Kdtree, please first plant the kdtree!"
        );
      }

      // good to go. ask the tree
      const int vote = 1;
      int id[vote];
      double distance[vote];
      kdtree_->knnSearch(goal,
                         vote /* closest neighbor */,
                         &id[0],
                         &distance[0]);

      using ReturnType = SplineType::ParametricCoordinate_;

      return std::move(
          grid_points_.IndexToParametricCoordinate<ReturnType>(id[0])
      );
    } // if initial_guess == ...
  }

  /*!
   * Builds RHS and fills "eyeders" matrix, which is also required in LHS.
   *
   * "eyeders" refer to collection of spline derivatives where only one
   * parametric dimension is set to one and everythin else zero, which results
   * in "eye"-like vectorized derivative query
   *
   * @params[in] guess current parametric coordinate guess
   * @params[in] difference result of `QueryMinusGuess()`
   * @params[out] eyeders
   * @params[out] rhs
   */
  void FillEyeDersAndRhs(
      const SplineType::ParametricCoordinate_& guess,
      const PArrayD_& difference,
      PxDMatrixD_& eyeders,
      PArrayD_& rhs) const {
    // btw, RHS is what's internally called as "df_dxi"
    double df_dxi_i; // temporary variable to hold sum
    int j;
    SplineType::Derivative_ derivative_query;

    for (int i{}; i < SplineType::kParaDim; ++i) {
      // eyeder query formulation
      derivative_query.fill(splinelib::Derivative{0});
      ++derivative_query[i];
      // derivative evaluation
      auto const &derivative = spline_(guess, derivative_query);

      j = 0;
      df_dxi_i = 0.;
      for (const auto& d : derivative) {
        // cast, since d may be a NamedType.
        const double d_value = static_cast<double>(d);
        df_dxi_i += difference[j] * d_value;
        eyeders[i][j] = d_value;
        ++j;
      }
      rhs[i] = -2. * df_dxi_i; // apply minus here already!
    }
  }

  /*!
   * Builds LHS
   *
   * @params[in] guess
   * @params[in] difference
   * @params[in] eyeders
   * @params[out] lhs
   */
  void FillLhs(
      const SplineType::ParametricCoordinate_& guess,
      const PArrayD_& difference,
      const PxDMatrixD_& eyeders,
      PxPMatrixD_& lhs) const {

    const PxPMatrixD_ eyedersAAt = AAt(eyeders);

    SplineType::Derivative_ derivative_query;
    for (int i{}; i < SplineType::kParaDim; ++i) {
      for (int j{i}; j < SplineType::kParaDim; ++j) {
        derivative_query.fill(splinelib::Derivative{0});
        ++derivative_query[i];
        ++derivative_query[j];

        auto const &derivative = spline_(guess, derivative_query);

        lhs[i][j] = 2 * (Dot(difference, derivative) + eyedersAAt[i][j]);

        if (i != j) lhs[j][i] = lhs[i][j]; // fill symmetric part
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
  SplineType::ParametricCoordinate_ FindNearestParametricCoordinate(
      const double* query,
      InitialGuess initial_guess,
      double tolerance = 1e-12) {

    PxPMatrixD_ lhs;
    PArrayD_ rhs;
    PArrayD_ delta_guess;
    DArrayD_ difference;
    PxDMatrixD_ eyeders;
    PArrayI_ clipped{}; /* clip status after most recent update */
    PArrayI_ previous_clipped{};
    PArrayI_ solver_skip_mask{}; /* tell solver to skip certain entry */
    bool solver_skip_mask_activated = false;

    // search_bounds is parametric bounds.
    const auto search_bounds =
        splinepy::splines::GetParametricBounds(splines_);

    SplineType::ParametricCoordinate_ current_guess =
        MakeInitialGuess(initial_guess, query);

    // Be optimistic and check if initial guess was awesome
    // If it is awesome, return and have feierabend
    QueryMinusGuess(query, current_guess, difference);
    if (splinepy::utils::NormL2(difference) < tolerance) {
      return std::move(current_guess);
    }

    const int max_iteration = SplineType::kParaDim * 20;
    double previous_norm{}, current_norm;

    // newton iteration
    for (int i{}; i < max_iteration; ++i) {
      // build systems to solve
      FillEyeDersAndRhs(current_guess, difference, eyeders, rhs);
      FillLhs(current_guess, difference, eyeders, lhs);

      // solve and update
      // 1. set solver skip mask if clipping happend twice at the same place.
      if (previous_clipped == clipped && !solver_skip_mask_activated) {
        solver_skip_mask = clipped;
        solver_skip_mask_activated = true;
        // if skip mask is on for all entries, return now.
        if (splinepy::utils::NonZeros(solver_skip_mask) ==
                SplineType::kParaDim) {
          // current_guess should be clipped at this point.
          return std::move(current_guess);
        }
      }

      // GaussWithPivot may swap and modify enties of all the input
      // -> can't use lhs and rhs afterwards, and we don't need them.
      // -> solver_skip_mask and delta_guess is reordered to rewind swaps
      splinepy::utils::GaussWithPivot(lhs, rhs, solver_skip_mask, delta_guess);

      // Update
      splinepy::utils::AddSecondToFirst(current_guess, delta_guess);
      // Clip
      previous_clipped = std::move(clipped);
      splinepy::utils::Clip(search_bounds, current_guess, clipped);
      // Converged?
      current_norm = NormL2(lhs);
      if (std::abs(previous_norm - current_norm) < tolernace) break;

      // prepare next round
      QueryMinusGuess(query, current_guess, difference);
      previous_norm = current_norm
    }
    return std::move(current_guess);
  }

protected:
  SplineType const &spline_;
  // kdtree related variables
  GridPoints_ grid_points_;
  PArrayI_ sampled_resolutions_;
  bool kdtree_planted_ = false;
  Coordinates_ coordinates_;
  Cloud_ coordinates_cloud_;
  Tree_ kdtree_;

};

} /* namespace proximity */
