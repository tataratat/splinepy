/*
MIT License

Copyright (c) 2021 Jaewook Lee

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "splinepy/proximity/proximity.hpp"
#include "splinepy/proximity/slsqp/slsqp.h"
#include "splinepy/splines/helpers/properties.hpp"
#include "splinepy/utils/nthreads.hpp"
#include "splinepy/utils/print.hpp"

namespace splinepy::proximity {

void SearchData::SLSQPData::Setup(SearchData& d) {
  // initialize options
  m = 0;
  meq = 0;
  la = 1;
  n = d.para_dim;
  x = d.current_guess.data(); //
  xl = d.lower_bound.data();
  xu = d.upper_bound.data();
  f = d.status.J; //
  c = nullptr;    // no constraints

  // allocate and set g
  g_space.Reallocate(n + 1);
  g_space[n] = 0.0;
  g = g_space.data(); // before calling, ccall CopyDerivativce() to fill this

  a = nullptr;
  // since SLSQP uses true J value, we adjust the tolerance accordingly
  acc = d.options.tolerance * d.options.tolerance * 0.5;
  iter = d.options.max_iter;
  mode = 0; // initialize

  // set workspace
  const int n1_ = n + 1;
  const int mineq = m - meq + n1_ + n1_;
  const int w_size = (3 * n1_ + m) * (n1_ + 1) + (n1_ - meq + 1) * (mineq + 2)
                     + 2 * mineq + (n1_ + mineq) * (n1_ - meq) + 2 * meq + n1_
                     + (((n + 1) * n) / 2) + 2 * m + 3 * n + 3 * n1_ + 1;

  // set w and lw
  w_space.Reallocate(w_size);
  w = w_space.data();
  l_w = w_size;

  // set jw and l_jw
  jw_space.Reallocate(mineq);
  jw = jw_space.data();
  l_jw = mineq;

  // internal states
  alpha = 0.0;
  f0 = 0.0;
  gs = 0.0;
  h1 = 0.0;
  h2 = 0.0;
  h3 = 0.0;
  h4 = 0.0;
  t = 0.0;
  t0 = 0.0;
  tol = 0.0;
  iexact = 0;
  incons = 0;
  ireset = 0;
  itermx = 0;
  line = 0;
  n1 = 0;
  n2 = 0;
  n3 = 0;
}

void SearchData::Setup(const int para_dim_,
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

void Proximity::PlantNewKdTree(const int* resolutions, const int n_thread) {
  const int para_dim = spline_.SplinepyParaDim();
  const int dim = spline_.SplinepyDim();

  // get parametric bounds
  RealArray parametric_bounds(para_dim * 2);
  spline_.SplinepyParametricBounds(parametric_bounds.data());

  // create grid_points
  grid_points_.SetUp(para_dim, parametric_bounds.data(), resolutions);

  const int n_queries = grid_points_.Size();

  // reallocate sampled spline
  sampled_spline_.Reallocate(n_queries * dim);

  // lambda function to allow n-thread execution
  auto sample_coordinates = [&](const int begin, const int end, int) {
    RealArray query(para_dim);
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

void Proximity::GuessMinusQuery(const RealArray& guess,
                                const ConstRealArray& query,
                                RealArray& difference) const {
  // evaluate guess and sett to difference
  spline_.SplinepyEvaluate(guess.data(), difference.data());
  // subtract query from evaluated guess
  difference.Subtract(query);
}

void Proximity::GuessMinusQuery(const RealArray& guess,
                                const ConstRealArray& query,
                                RealArray& difference,
                                RealArray& guess_phys) const {
  spline_.SplinepyEvaluate(guess.data(), guess_phys.data());
  splinepy::utils::Subtract(guess_phys, query, difference);
}

void Proximity::MakeInitialGuess(const ConstRealArray& goal,
                                 RealArray& guess) const {
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

double Proximity::J(const RealArray& guess,
                    const ConstRealArray& query,
                    RealArray& difference,
                    RealArray& guess_phys) const {
  GuessMinusQuery(guess, query, difference, guess_phys);
  return 0.5 * difference.NormL2Squared();
}

void Proximity::dJdu(const RealArray& guess,
                     const RealArray& difference,
                     RealArray& djdu,
                     RealArray2D& spline_gradient) const {
  const int para_dim = guess.size();
  const int dim = difference.size();

  IndexArray derivative_query(para_dim);
  derivative_query.Fill(0);

  // this is just to view and apply inner product
  RealArray gradient_row_view;
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

void Proximity::d2Jdu2(const RealArray& guess,
                       const RealArray& difference,
                       const RealArray2D& spline_gradient_AAt,
                       RealArray2D& d2jdu2,
                       RealArray3D& spline_hessian) const {
  const int para_dim = guess.size();
  const int dim = difference.size();

  IndexArray derivative_query(para_dim);
  derivative_query.Fill(0);

  // derivative result is a view to the hessian array
  RealArray derivative;
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
                      data.difference,
                      data.current_phys);
  }

  if (depth > 0) {
    dJdu(data.current_guess, data.difference, data.djdu, data.spline_gradient);
  }

  if (depth > 1) {
    data.spline_gradient.AAt(data.spline_gradient_AAt);
    d2Jdu2(data.current_guess,
           data.difference,
           data.spline_gradient_AAt,
           data.d2jdu2,
           data.spline_hessian);
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
  spline_.SplinepyParametricBounds(data.search_bounds.data());

  // set search bounds to para bounds and early exit
  if (!tight) {
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
  auto& lhs = data.Lhs(data.d2jdu2);
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
    UpdateAndClip(data);

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
  auto& lhs = data.Lhs(data.spline_gradient_AAt);
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
  data.LM.Setup(data);
  while (i < data.options.max_iter) {
    if (data.IsConverged()) {
      break;
    }
    lhs.ResetRowOrder();
    lhs.Solve(rhs, data.delta_guess);

    // update the save original_guess before updating
    data.LM.prev_guess = data.current_guess;
    UpdateAndClip(data);

    // get actual delta guess, in case it is clipped
    if (data.clipped.NonZeros()) {
      data.delta_guess = data.LM.prev_guess;
      data.delta_guess.Subtract(data.current_guess);
    }

    // compute two different metrics and take bigger one
    auto& metric = data.LM.Metric(data.spline_gradient, data.delta_guess);
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
        // if lambda is between 2^13 and 2^14, we exit
        if (data.LM.lambda > 1.e5) {
          // status needs re-evaluation at updated location
          ComputeStatus(data, true);
          break;
        }

        // we continue -> set current_guess back to previous guess
        data.current_guess = data.LM.prev_guess;
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

void Proximity::PrepareIterationSlsqp(SearchData& d) const {
  // precautious step based on
  // https://github.com/scipy/scipy/issues/11403
  d.current_guess.Clip(d.lower_bound, d.upper_bound, d.clipped);

  // mode 0: this is first
  if (d.SLSQP.mode == 0) {
    ComputeCostAndDerivatives(d, 1);
    ComputeStatus(d, false); // strictly not required
    d.SLSQP.CopyDerivative(d.djdu);
    d.SLSQP.f = d.status.J;
  } else if (d.SLSQP.mode == 1) {
    ComputeCostAndDerivatives(d, 0);
    d.SLSQP.f = d.status.J;
  } else if (d.SLSQP.mode == -1) {
    ComputeCostAndDerivatives(d, 1);
    d.SLSQP.CopyDerivative(d.djdu);
  }
}

void Proximity::Slsqp(SearchData& d) const {
  d.status.stop_iteration = 0;

  d.SLSQP.Setup(d);
  PrepareIterationSlsqp(d);

  ComputeStatus(d, false);

  // set convergence goal based on the first round
  d.options.convergence_goal = ConvergenceGoal(d.status.convergence_norm,
                                               d.options.tolerance,
                                               d.options.tolerance);

  int i{};
  while (true) {

    // call slsqp
    slsqp(&d.SLSQP.m,
          &d.SLSQP.meq,
          &d.SLSQP.la,
          &d.SLSQP.n,
          d.SLSQP.x,
          d.SLSQP.xl,
          d.SLSQP.xu,
          &d.SLSQP.f,
          d.SLSQP.c,
          d.SLSQP.g,
          d.SLSQP.a,
          &d.SLSQP.acc,
          &d.SLSQP.iter,
          &d.SLSQP.mode,
          d.SLSQP.w,
          &d.SLSQP.l_w,
          d.SLSQP.jw,
          &d.SLSQP.l_jw,
          &d.SLSQP.alpha,
          &d.SLSQP.f0,
          &d.SLSQP.gs,
          &d.SLSQP.h1,
          &d.SLSQP.h2,
          &d.SLSQP.h3,
          &d.SLSQP.h4,
          &d.SLSQP.t,
          &d.SLSQP.t0,
          &d.SLSQP.tol,
          &d.SLSQP.iexact,
          &d.SLSQP.incons,
          &d.SLSQP.ireset,
          &d.SLSQP.itermx,
          &d.SLSQP.line,
          &d.SLSQP.n1,
          &d.SLSQP.n2,
          &d.SLSQP.n3);

    PrepareIterationSlsqp(d);

    if (std::abs(d.SLSQP.mode) != 1) {
      break;
    }

    d.status.stop_iteration = ++i;
  }

  // compute status for further pipeline
  // once slsqp finds solution within the tolernace, it evaluates g once more
  // so we probably don't need to recompute.
  // just to be safe.
  ComputeStatus(d, true);
}

void Proximity::VerboseQuery(
    const double* query,
    const double& tolerance,
    const int& max_iterations,
    const bool tight_bounds,
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

  // we want to keep the best result
  ComputeCostAndDerivatives(data, 0); // computes J
  double best_J = data.status.J;
  RealArray best_guess(data.initial_guess);
  auto keep_best_and_did_it_improve = [&best_J, &best_guess, &data]() -> bool {
    // current guess is better - use the new one
    if (best_J > data.status.J) {
      best_guess = data.current_guess;
      best_J = data.status.J;
      return true;
    }

    // previous guess was better - set back
    data.current_guess = best_guess;
    best_J = data.status.J;
    return false;
  };

  // return preps
  auto fill_hessian_and_metrics =
      [&data, &distance, &convergence_norm, this](const bool compute) {
        if (!data.spline_hessian.OwnsData()) {
          if (compute) {
            FillHessian(spline_, data.current_guess, data.spline_hessian);
          } else {
            splinepy::utils::CopyUpperToLowerTriangle(data.spline_hessian);
          }
        }

        // metrics
        distance = data.status.distance;
        convergence_norm = data.status.convergence_norm;
      };

  // set bounds
  FindSearchBound(data, tight_bounds);

  // we first start with Newton
  Newton(data);

  // newton converged - just fill lower triangle of the second_derivative
  if (data.IsConverged()) {
    // after newton, we just need to copy upper to lower triangle
    fill_hessian_and_metrics(false);
    return;
  }

  // newton didn't converge. prepare SLSQP
  keep_best_and_did_it_improve();
  data.options.max_iter *= 5;
  Slsqp(data);
  // check
  if (data.IsConverged() || data.SLSQP.mode == 0) {
    fill_hessian_and_metrics(true);
    return;
  }

  // Newton and Slsqp didn't work. Try LevenbergMarquart
  keep_best_and_did_it_improve();
  // set max iteration higher
  data.options.max_iter *= 5;
  LevenbergMarquart(data);

  // if it is improved, recompute all
  if (keep_best_and_did_it_improve()) {
    ComputeStatus(data, true); // this should fill grad
  }

  // finally, hessian and metrics
  fill_hessian_and_metrics(true);
}

void FillHessian(const splinepy::splines::SplinepyBase& spline,
                 const RealArray& at,
                 RealArray3D& hess) {
  const int para_dim = hess.Shape()[1];
  const int dim = hess.Shape()[2];

  IndexArray derivative_query(para_dim);
  derivative_query.Fill(0);

  // derivative result is a view to the hessian array
  RealArray derivative;
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
