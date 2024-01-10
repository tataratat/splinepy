#include "splinepy/proximity/proximity.hpp"
#include "splinepy/splines/helpers/properties.hpp"
#include "splinepy/utils/nthreads.hpp"
#include "splinepy/utils/print.hpp"

namespace splinepy::proximity {

void Proximity::PlantNewKdTree(const int* resolutions, const int n_thread) {
  const int para_dim = spline_.SplinepyParaDim();
  const int dim = spline_.SplinepyDim();

  // get parametric bounds
  RealArray_ parametric_bounds(para_dim * 2);
  spline_.SplinepyParametricBounds(parametric_bounds.data());

  // create grid_points
  grid_points_.SetUp(para_dim, parametric_bounds.data(), resolutions);

  const int n_queries = grid_points_.Size();

  // reallocate sampled spline
  sampled_spline_.Reallocate(n_queries * dim);

  // lambda function to allow n-thread execution
  auto sample_coordinates = [&](const int begin, const int end, int) {
    RealArray_ query(para_dim);
    double* query_data = query.data();

    for (int i{begin}; i < end; ++i) {
      grid_points_.IdToGridPoint(i, query_data);
      spline_.SplinepyEvaluate(query_data, &sampled_spline_[i * dim]);
    }
  };

  // n-thread execution for sampling
  splinepy::utils::NThreadExecution(sample_coordinates, n_queries, n_thread);

  // nanoflann supports concurrent build
  nanoflann::KDTreeSingleIndexAdaptorParams params{};
  params.n_thread_build = static_cast<
      decltype(nanoflann::KDTreeSingleIndexAdaptorParams::n_thread_build)>(
      (n_thread < 0) ? 0 : n_thread);

  // plant a new tree
  cloud_ = std::make_unique<Cloud_>(sampled_spline_.data(),
                                    sampled_spline_.size(),
                                    dim);
  kdtree_ = std::make_unique<Tree_>(dim, *cloud_, params);
}

void Proximity::GuessMinusQuery(const RealArray_& guess,
                                const ConstRealArray_& query,
                                RealArray_& difference) const {
  // evaluate guess and sett to difference
  spline_.SplinepyEvaluate(guess.data(), difference.data());
  // subtract query from evaluated guess
  difference.Subtract(query);
}

void Proximity::GuessMinusQuery(const RealArray_& guess,
                                const ConstRealArray_& query,
                                RealArray_& guess_phys,
                                RealArray_& difference) const {
  spline_.SplinepyEvaluate(guess.data(), guess_phys.data());
  splinepy::utils::Subtract(guess_phys, query, difference);
}

void Proximity::MakeInitialGuess(const ConstRealArray_& goal,
                                 RealArray_& guess) const {
  if (!kdtree_) {
    // hate to be aggressive, but here it is.
    splinepy::utils::PrintAndThrowError("to use InitialGuess::Kdtree option,"
                                        "please first plant a kdtree.");
  }

  // good to go. ask the tree
  int id;
  double distance;
  kdtree_->knnSearch(goal.data(), 1 /* closest neighbor */, &id, &distance);

  grid_points_.IdToGridPoint(id, guess.data());
}

void Proximity::FillSplineGradientAndRhs(const RealArray_& guess,
                                         const RealArray_& difference,
                                         RealArray2D_& spline_gradient,
                                         RealArray_& rhs) const {
  const int para_dim = guess.size();
  const int dim = difference.size();

  IndexArray_ derivative_query(para_dim);
  derivative_query.Fill(0);

  // this is just to view and apply inner product
  RealArray_ gradient_row_view;
  gradient_row_view.SetShape(dim);

  for (int i{}; i < para_dim; ++i) {
    // set query - this should be all zero here
    derivative_query[i] = 1;

    // set row view
    gradient_row_view.SetData(&spline_gradient(i, 0));

    // derivative evaluation
    spline_.SplinepyDerivative(guess.data(),
                               derivative_query.data(),
                               gradient_row_view.data());

    // fill rhs_i and apply minus here already!
    rhs[i] = -2. * difference.InnerProduct(gradient_row_view);

    // reset query to zero
    derivative_query[i] = 0;
  }
}

void Proximity::FillLhs(const RealArray_& guess,
                        const RealArray_& difference,
                        const RealArray2D_& spline_gradient_AAt,
                        RealArray3D_& spline_hessian,
                        RealArray2D_& lhs) const {
  const int para_dim = guess.size();
  const int dim = difference.size();

  IndexArray_ derivative_query(para_dim);
  derivative_query.Fill(0);

  // derivative result is a view to the hessian array
  RealArray_ derivative;
  derivative.SetShape(dim);

  // lambda to compute each element
  auto compute = [&](const int& i, const int& j) {
    // adjust derivative query
    ++derivative_query[i];
    ++derivative_query[j];

    // adjust view on hessian array
    derivative.SetData(&spline_hessian(i, j, 0));

    // compute
    spline_.SplinepyDerivative(guess.data(),
                               derivative_query.data(),
                               derivative.data());

    // fill lhs_ij
    lhs(i, j) =
        2. * (difference.InnerProduct(derivative) + spline_gradient_AAt(i, j));

    // reset derivative query
    --derivative_query[i];
    --derivative_query[j];
  };

  for (int i{}; i < para_dim; ++i) {
    for (int j{i + 1}; j < para_dim; ++j) {
      // upper triangle without diagonal
      compute(i, j);
      // copy to the lower
      lhs(j, i) = lhs(i, j);
    }
    // diagonal
    compute(i, i);
  }
}

void Proximity::VerboseQuery(
    const double* query,
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

  const int para_dim = spline_.SplinepyParaDim();
  const int dim = spline_.SplinepyDim();

  // view arrays - we will use this memory for IO
  // this avoid unnecessary copy, but if it slows down the process
  // significantly we will just alloc and copy at the end
  // RealArray_ phys_query(query, dim);
  ConstRealArray_ phys_query(query, dim);
  RealArray_ current_guess(final_guess, para_dim);
  RealArray_ current_phys(nearest, dim);
  RealArray_ difference(nearest_minus_query, dim);
  RealArray2D_ spline_gradient(first_derivatives, para_dim, dim);
  RealArray3D_ spline_hessian(second_derivatives, para_dim, para_dim, dim);

  // allocate aux real arrays
  RealArray2D_ lhs(para_dim, para_dim);
  SystemMatrix system(lhs);
  RealArray_ rhs(para_dim);
  RealArray_ delta_guess(para_dim);
  RealArray2D_ spline_gradient_AAt(para_dim, para_dim);
  RealArray2D_ search_bounds(2, para_dim);

  // get pointers to beginning of each bound
  RealArray_ lower_bound(search_bounds.begin(), para_dim);
  RealArray_ upper_bound(search_bounds.begin() + para_dim, para_dim);

  // allocate index arrays
  IndexArray_ clipped(para_dim);
  IndexArray_ previous_clipped(para_dim);

  // search_bounds is parametric bounds here
  spline_.SplinepyParametricBounds(search_bounds.data());

  // initial guess
  MakeInitialGuess(phys_query, current_guess);

  // Let's try aggressive search bounds
  if (aggressive_bounds) {
    // you need to be sure that you have sampled your spline fine enough
    for (int i{}; i < para_dim; ++i) {
      // adjust lower (0) and upper (1) bounds aggressively
      // but of course, not so aggressive that it is out of bound.
      search_bounds(0, i) =
          std::max(search_bounds(0, i),
                   current_guess[i] - grid_points_.step_size_[i]);
      search_bounds(1, i) =
          std::min(search_bounds(1, i),
                   current_guess[i] + grid_points_.step_size_[i]);
    }
  }

  // get initial status - call both rhs and lhs to fill gradients
  GuessMinusQuery(current_guess, phys_query, current_phys, difference);
  FillSplineGradientAndRhs(current_guess, difference, spline_gradient, rhs);
  spline_gradient.AAt(spline_gradient_AAt);
  FillLhs(current_guess, difference, spline_gradient_AAt, spline_hessian, lhs);

  distance = difference.NormL2();
  convergence_norm = rhs.NormL2();

  // for now, we take same value for rel_tol
  const double norm_goal = std::max(convergence_norm * tolerance, tolerance);
  // get maxiteration
  const int max_iter = max_iterations < 0 ? para_dim * 20 : max_iterations;
  // newton iterations
  for (int i{}; i < max_iter; ++i) {
    // norm and distance check for convergence
    if (convergence_norm < norm_goal || distance < tolerance) {
      break;
    }

    // solve systems using gauss elimination with partial pivoting
    // this will alter lhs in place, but shouldn't reorder rows inplace
    system.Solve(rhs, delta_guess);

    // currently we just clip at the bounds
    // we need a more sophisticated update for stability
    current_guess.Add(delta_guess);
    current_guess.Clip(lower_bound, upper_bound, clipped);

    // evaluate cost at current guess
    GuessMinusQuery(current_guess, phys_query, current_phys, difference);
    distance = difference.NormL2();

    // assemble for next iteration
    // this also count as rhs/lhs at current guess, which fills derivatives
    // for return
    FillSplineGradientAndRhs(current_guess, difference, spline_gradient, rhs);
    spline_gradient.AAt(spline_gradient_AAt);
    FillLhs(current_guess,
            difference,
            spline_gradient_AAt,
            spline_hessian,
            lhs);
    convergence_norm = rhs.NormL2();
  }

  // before returning, we just need to fill symmetric part of hessian
  for (int i{}; i < spline_hessian.Shape()[0]; ++i) {
    for (int j{i + 1}; j < spline_hessian.Shape()[1]; ++j) {
      std::copy_n(&spline_hessian(i, j, 0), dim, &spline_hessian(j, i, 0));
    }
  }
}

} // namespace splinepy::proximity
