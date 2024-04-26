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

void Proximity::MakeInitialGuess(SearchData& data) const {
  MakeInitialGuess(data.phys_query, data.initial_guess);
  data.current_guess = data.initial_guess;
}

void Proximity::UpdateAndClip(SearchData& data) const {
  data.current_guess.Add(data.delta_guess);
  data.current_guess.Clip(data.lower_bound, data.upper_bound, data.clipped);
};

double Proximity::J(const RealArray_& guess,
                    const ConstRealArray_& query,
                    RealArray_& guess_phys,
                    RealArray_& difference) const {
  GuessMinusQuery(guess, query, guess_phys, difference);
  return 0.5 * difference.NormL2Squared();
}

void Proximity::dJdu(const RealArray_& guess,
                     const RealArray_& difference,
                     RealArray2D_& spline_gradient,
                     RealArray_& djdu) const {
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

    // note that there's no negative here, thus not rhs yet
    djdu[i] = difference.InnerProduct(gradient_row_view);

    // reset query to zero
    derivative_query[i] = 0;
  }
}

void Proximity::d2Jdu2(const RealArray_& guess,
                       const RealArray_& difference,
                       const RealArray2D_& spline_gradient_AAt,
                       RealArray3D_& spline_hessian,
                       RealArray2D_& d2jdu2) const {
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
    d2jdu2(i, j) =
        difference.InnerProduct(derivative) + spline_gradient_AAt(i, j);

    // reset derivative query
    --derivative_query[i];
    --derivative_query[j];
  };

  for (int i{}; i < para_dim; ++i) {
    for (int j{i + 1}; j < para_dim; ++j) {
      // upper triangle without diagonal
      compute(i, j);
      // copy to the lower
      d2jdu2(j, i) = d2jdu2(i, j);
    }
    // diagonal
    compute(i, i);
  }
}

void Proximity::ComputeCostAndDerivatives(SearchData& data, int depth) const {
  // magic negative will set to the highest depth we use.
  depth = (depth < 0) ? 2 : depth;

  if (depth > -1) {
    data.status.J = J(data.current_guess,
                      data.phys_query,
                      data.current_phys,
                      data.difference);
  }

  if (depth > 0) {
    dJdu(data.current_guess, data.difference, data.spline_gradient, data.djdu);
  }

  if (depth > 1) {
    data.spline_gradient.AAt(data.spline_gradient_AAt);
    d2Jdu2(data.current_guess,
           data.difference,
           data.spline_gradient_AAt,
           data.spline_hessian,
           data.d2jdu2);
  }
}

void Proximity::ComputeStatus(SearchData& data,
                              const bool recompute_costs) const {
  if (recompute_costs) {
    ComputeCostAndDerivatives(data, 1);
  }

  // fill status
  data.status.squared_distance = data.difference.NormL2Squared();
  data.status.distance = std::sqrt(data.status.squared_distance);
  data.status.J = 0.5 * data.status.squared_distance;
  data.status.convergence_norm = data.djdu.NormL2();
}

void Proximity::FindSearchBound(SearchData& data, const bool tight) const {
  if (!tight) {
    // set search bounds to para bounds and early exit
    spline_.SplinepyParametricBounds(data.search_bounds.data());
    return;
  }

  // you need to be sure that you have sampled your spline fine enough
  for (int i{}; i < data.para_dim; ++i) {
    // adjust lower (0) and upper (1) bounds aggressively
    // but of course, not so aggressive that it is out of bound.
    data.search_bounds(0, i) =
        std::max(data.search_bounds(0, i),
                 data.current_guess[i] - grid_points_.step_size_[i]);
    data.search_bounds(1, i) =
        std::min(data.search_bounds(1, i),
                 data.current_guess[i] + grid_points_.step_size_[i]);
  }
}

double Proximity::LineSearchUpdate(SearchData& data) const {
  RealArray_ orig_guess = data.current_guess;
  RealArray_ orig_delta = data.delta_guess;

  const double q1 = data.status.J; // this should be from previous iteration

  // Try current full step
  UpdateAndClip(data);
  // get actual delta guess, in case it is clipped
  // we will consider this as our full step
  if (data.clipped.NonZeros()) {
    splinepy::utils::Subtract(orig_guess, data.current_guess, data.delta_guess);
  }

  // get residual of full step
  ComputeCostAndDerivatives(data, 1);
  ComputeStatus(data, false);
  const double q3 = data.status.J;

  // get one from half step and get residual
  data.current_guess.Add(-.5, data.delta_guess);
  ComputeCostAndDerivatives(data, 1);
  ComputeStatus(data, false);

  const double q2 = data.status.J;

  const double eps = (3.0 * q1 - 4.0 * q2 + q3) / (4.0 * (q1 - 2.0 * q2 + q3));

  double beta;
  if ((q1 - 2.0 * q2 + q3) > 0 && eps > 0 && eps < 1) {
    beta = eps;
  } else if (q3 < q1) {
    beta = 1.0;
  } else {
    beta = 0.05;
  }

  // use beta to update
  splinepy::utils::Add(orig_guess, beta, data.delta_guess, data.current_guess);

  return beta;
}

void Proximity::PrepareIterationNewton(SearchData& data) const {

  // newton requires hessian
  ComputeCostAndDerivatives(data, 2);

  // to be ready for solve, flip sign for rhs
  data.djdu.FlipSign();
}

void Proximity::Newton(SearchData& data) const {
  // initialize stop_iteration
  data.status.stop_iteration = 0;

  /// set lhs matrix to d2jdu2
  auto& lhs = data.LhsMatrix(data.d2jdu2);
  auto& rhs = data.djdu;

  // first round
  PrepareIterationNewton(data);
  ComputeStatus(data, false); // compute status based on current values

  // set convergence goal based on the first round
  data.options.convergence_goal = ConvergenceGoal(data.status.convergence_norm,
                                                  data.options.tolerance,
                                                  data.options.tolerance);

  int i{};
  while (i < data.options.max_iter) {
    if (data.IsConverged()) {
      break;
    }
    // solve systems using gauss elimination with partial pivoting
    // this will alter lhs in place, but shouldn't reorder rows inplace
    lhs.ResetRowOrder();
    lhs.Solve(rhs, data.delta_guess);

    // x_n+1 = x_n + delta_x
    // UpdateAndClip(data);
    LineSearchUpdate(data);

    // Prepare next step and compute status
    PrepareIterationNewton(data);
    ComputeStatus(data, false);

    data.status.stop_iteration = ++i;
  }
}

void Proximity::PrepareIterationLevenbergMarquart(SearchData& data) const {
  // LM requires Jac (pronounced similar to Jacques)
  ComputeCostAndDerivatives(data, 1);

  // prepare rhs
  data.djdu.FlipSign();

  // prepare lhs
  data.spline_gradient.AAt(data.spline_gradient_AAt);

  // modify lhs based on current value in data.LM
  double diag_contribution{data.LM.lambda};
  const bool modified = data.LM.modified_update;
  if (modified) {
    diag_contribution += 1.0;
  }

  for (int i{}; i < data.para_dim; ++i) {
    if (modified) {
      data.spline_gradient_AAt(i, i) *= diag_contribution;
    } else {
      data.spline_gradient_AAt(i, i) += diag_contribution;
    }
  }
}

void Proximity::LevenbergMarquart(SearchData& data) const {
  // initialize stop_iteration
  data.status.stop_iteration = 0;

  // for lhs we will use AAt matrix. Modifications may be made to diagonal
  // entries
  auto& lhs = data.LhsMatrix(data.spline_gradient_AAt);
  auto& rhs = data.djdu;

  // before first round, initialize some options
  data.LM.lambda = 1.;
  data.LM.lower_bound = .25;
  data.LM.upper_bound = .75;

  PrepareIterationLevenbergMarquart(data);
  ComputeStatus(data, false);

  // set initial previous value
  data.LM.prev_J = data.status.J;

  // set convergence goal based on the first round
  data.options.convergence_goal = ConvergenceGoal(data.status.convergence_norm,
                                                  data.options.tolerance,
                                                  data.options.tolerance);

  int i{};
  RealArray_ prev_guess(data.para_dim);
  RealArray_ metric(data.para_dim);
  while (i < data.options.max_iter) {
    if (data.IsConverged()) {
      break;
    }
    lhs.ResetRowOrder();
    lhs.Solve(rhs, data.delta_guess);

    // update the save original_guess before updating
    prev_guess = data.current_guess;
    UpdateAndClip(data);

    // get actual delta guess, in case it is clipped
    if (data.clipped.NonZeros()) {
      data.delta_guess = prev_guess;
      data.delta_guess.Subtract(data.current_guess);
    }

    // compute two different metrics and take bigger one
    data.LM.Metric(data.spline_gradient, data.delta_guess, metric);
    const double metric2 = metric.InnerProduct(data.difference);

    metric.Add(data.difference);
    const double metric_norm = .5 * metric.NormL2Squared();

    // compute || F(x + dX) ||^2 then rho
    ComputeCostAndDerivatives(data, 0); // computes J and everything on the way
    const double rho =
        (data.LM.prev_J - data.status.J) / (data.LM.prev_J - metric_norm);
    const double rho2 = -(data.LM.prev_J - data.status.J) / metric2;

    const double chosen_rho = std::max(rho, rho2);

    // If the metric is smaller than the lower bound, the update is
    // rejected and the penelization is increased
    if (chosen_rho < data.LM.lower_bound) {
      // This is @jzwar's practical trick:
      // Increase lambda if current value is worse than previous
      if (data.status.J > data.LM.prev_J) {
        data.LM.lambda *= 2.0;
        // if it gets too big, exit
        // TODO: better "too big" number
        if (data.LM.lambda > 1.e5) {
          break;
        }

        // we continue -> set current_guess back to previous guess
        data.current_guess = prev_guess;
        PrepareIterationLevenbergMarquart(data);
        ComputeStatus(data, false);
        continue;
      }
    }

    // if rho is too big, decrease lambda
    if (chosen_rho > data.LM.upper_bound) {
      data.LM.lambda *= 0.5;
    }

    // we are here, because lambda is not modified
    // take this value - current_guess stays the same,
    data.LM.prev_J = data.status.J;

    PrepareIterationLevenbergMarquart(data);
    ComputeStatus(data, false);

    data.status.stop_iteration = ++i;
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

  // get dims
  const int para_dim = spline_.SplinepyParaDim();
  const int dim = spline_.SplinepyDim();

  // initialize aux data for the search
  SearchData data(para_dim,
                  dim,
                  tolerance,
                  max_iterations,
                  query,
                  final_guess,
                  nearest,
                  nearest_minus_query,
                  first_derivatives,
                  second_derivatives);

  // initial guess and set it as current guess
  MakeInitialGuess(data);

  // set bounds
  FindSearchBound(data, aggressive_bounds);

  // we first start with Newton
  Newton(data);

  // newton converged - just fill lower triangle of the second_derivative
  if (data.IsConverged()) {
    // we only need to fill data if hessian array belongs to the function caller
    // otherwise, it means it's created within this function and so don't care.
    if (!data.spline_hessian.OwnsData()) {
      splinepy::utils::CopyUpperToLowerTriangle(data.spline_hessian);
    }
    return;
  }

  // // Newton didn't work. Try LevenbergMarquart
  // // reset current_guess to initial_guess as it expects that
  // data.current_guess = data.initial_guess;
  // LevenbergMarquart(data);
  // FillHessian(spline_, data.current_guess, data.spline_hessian);
}

void FillHessian(const splinepy::splines::SplinepyBase& spline,
                 const Proximity::RealArray_& at,
                 Proximity::RealArray3D_& hess) {
  const int para_dim = hess.Shape()[1];
  const int dim = hess.Shape()[2];

  Proximity::IndexArray_ derivative_query(para_dim);
  derivative_query.Fill(0);

  // derivative result is a view to the hessian array
  Proximity::RealArray_ derivative;
  derivative.SetShape(dim);

  auto compute = [&](const int& i, const int& j) {
    // adjust derivative query
    ++derivative_query[i];
    ++derivative_query[j];

    // adjust view on hessian array
    derivative.SetData(&hess(i, j, 0));

    // compute
    spline.SplinepyDerivative(at.data(),
                              derivative_query.data(),
                              derivative.data());

    // reset derivative query
    --derivative_query[i];
    --derivative_query[j];
  };

  for (int i{}; i < para_dim; ++i) {
    for (int j{i + 1}; j < para_dim; ++j) {
      // upper triangle without diagonal
      compute(i, j);
      // copy derivative to lower
      std::copy_n(derivative.data(), dim, &hess(j, i, 0));
    }
    // diagonal
    compute(i, i);
  }
}

} // namespace splinepy::proximity
