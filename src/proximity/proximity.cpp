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

void Proximity::FillLhsLM(const RealArray_& guess,
                          const RealArray_& difference,
                          const RealArray2D_& spline_gradient_AAt,
                          RealArray2D_& lhs,
                          const double& lambda,
                          const bool& modified_marquart) const {
  const int para_dim = guess.size();

  for (int i{}; i < para_dim; ++i) {
    for (int j{i + 1}; j < para_dim; ++j) {
      // upper triangle without diagonal ( two is to match RHS from Newton)
      lhs(i, j) = 2 * spline_gradient_AAt(i, j);
      // copy to the lower
      lhs(j, i) = lhs(i, j);
    }
    // diagonal
    if (modified_marquart) {
      lhs(i, i) = 2 * spline_gradient_AAt(i, i) * (1 + lambda);
    } else {
      lhs(i, i) = 2 * (spline_gradient_AAt(i, i) + lambda);
    }
  }
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
  auto& lhs = data.LhsMatrix(data.spline_gradient_AAt);
  auto& rhs = data.djdu;

  // before first round, initialize some options
  data.LM.lambda = 1.;
  data.LM.lower_bound = .25;
  data.LM.upper_bound = .75;

  PrepareIterationLevenbergMarquart(data);
  ComputeStatus(data, false);

  data.LM.prev_J = data.status.J;

  int i{};
  RealArray_ prev_guess(data.para_dim);
  RealArray_ metric(data.para_dim);
  double alt_metric;
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
      // i think this should have "continue" here - double check!
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

  // Newton didn't work. Try LevenbergMarquart
  // reset current_guess to initial_guess as it expects that
  data.current_guess = data.initial_guess;
  LevenbergMarquart(data);

  FillHessian(spline_, data.current_guess, data.spline_hessian);

  // // Check if converged
  // if (convergence_norm > norm_goal && distance > tolerance) {
  //   // Did not converge, try LM
  //   std::cout << "START :";
  //   // Set back to starting point
  //   current_guess = initial_guess;
  //   GuessMinusQuery(current_guess, phys_query, current_phys, difference);

  //   // Set auxiliary values
  //   double lambda = 1.;
  //   double former_squared_distance = difference.NormL2Squared();
  //   RealArray_ current_difference = difference;

  //   // Coefficients for Levenberg-Marquart algorithm (empirical)
  //   const double lower_bound_LM = 0.25;
  //   const double upper_bound_LM = 0.75;

  //   // Calculate RHS before first iteration (remains unchanged by update)
  //   FillSplineGradientAndRhs(current_guess, difference, spline_gradient,
  //   rhs);

  //   for (int i_iteration{}; i_iteration < max_iter; ++i_iteration) {
  //     // norm and distance check for convergence
  //     if (convergence_norm < norm_goal || distance < tolerance) {
  //       std::cout << "SUCCESS" << std::endl;
  //       break;
  //     }
  //     // std::cout << "---------------\nIteration : " << i << std::endl;

  //     // Compute J^T J (Jacobian is transposed gradient)
  //     spline_gradient.AAt(spline_gradient_AAt);
  //     FillLhsLM(current_guess,
  //               current_difference,
  //               spline_gradient_AAt,
  //               lhs,
  //               lambda,
  //               false); // test true here

  //     // solve systems using gauss elimination with partial pivoting
  //     // this will alter lhs in place, but shouldn't reorder rows inplace
  //     system.Solve(rhs, delta_guess);

  //     // Decide whether or not to accept the next value by computing the
  //     metric
  //     // for rho
  //     RealArray_ next_guess = current_guess;
  //     next_guess.Add(delta_guess);
  //     next_guess.Clip(lower_bound, upper_bound, clipped);

  //     // Update delete to clipped delta
  //     delta_guess = next_guess;
  //     delta_guess.Subtract(current_guess);

  //     // Compute || F(x) + J^T * dX ||^2
  //     RealArray_ metric(delta_guess.size());

  //     // Inner product not supported for matrices
  //     // std::cout << "Delta Guess : ";
  //     // delta_guess.Print();
  //     // std::cout << "J^T : ";
  //     // spline_gradient.Print();
  //     // std::cout << " J^T dX :";
  //     for (int i{}; i < para_dim; ++i) {
  //       metric[i] = spline_gradient(i, 0) * delta_guess(0);
  //       for (int j{1}; j < dim; ++j) {
  //         metric[i] += spline_gradient(i, j) * delta_guess(j);
  //       }
  //     }
  //     // Konstantins metric
  //     const double alt_metric = metric.InnerProduct(current_difference);

  //     // metric.Print();
  //     metric.Add(current_difference);
  //     const double metric_norm = metric.NormL2Squared();

  //     // Compute || F(x + dX) ||^2
  //     GuessMinusQuery(next_guess, phys_query, current_phys, difference);
  //     const double distance_squared = difference.NormL2Squared();

  //     // Compute Rho
  //     const double rho = (former_squared_distance - distance_squared)
  //                        / (former_squared_distance - metric_norm);

  //     // Konstantin's notes
  //     const double rho_alt =
  //         -(former_squared_distance - distance_squared) / alt_metric;

  //     const double my_rho = std::max(rho, rho_alt);

  //     // // Debugging output
  //     std::cout << std::setw(15) << my_rho << " " << std::setw(15)
  //               << distance_squared << " " << std::setw(15) << metric_norm
  //               << std::setw(15) << former_squared_distance << " "
  //               << std::setw(15) << delta_guess.NormL2Squared() << " "
  //               << std::setw(15) << lambda << "  " << std::flush;

  //     // If the metric is smaller than the lower bound, the update is
  //     // rejected and the penelization is increased
  //     if ((my_rho < lower_bound_LM)
  //         // This part is not part of any notes, but it works better somehow.
  //         // KKs notes say nowhere that lambda should remain unchanged even
  //         if
  //         // the value is better than the previous. I added that here
  //         && (distance_squared > former_squared_distance)) {
  //       lambda *= 2.0;
  //       // std::cout << "+";
  //       if (lambda > 1e5) {
  //         break;
  //       }
  //       // std::cout << "Rejected ! New Lambda = " << lambda << std::endl;
  //       if (distance_squared > former_squared_distance) {
  //         std::cout << "Rejected" << std::endl;
  //         continue;
  //       }
  //     }

  //     std::cout << "Accepted" << std::endl;
  //     // If the metric is bigger than the upper bound the penelization is
  //     // divided by two
  //     if (my_rho > upper_bound_LM) {
  //       lambda *= 0.5;
  //       // std::cout << "-";
  //     } else {
  //       // std::cout << "=";
  //     }
  //     // std::cout << "Accepted ! New Lambda = " << lambda << std::endl;

  //     // If lower < rho < upper it is accepted and lambda is not modified
  //     // update values
  //     former_squared_distance = distance_squared;
  //     current_difference = difference;
  //     current_guess = next_guess;

  //     // evaluate cost at current guess
  //     distance = std::sqrt(former_squared_distance);

  //     // assemble for next iteration
  //     // this also count as rhs/lhs at current guess, which fills derivatives
  //     // for return
  //     FillSplineGradientAndRhs(current_guess, difference, spline_gradient,
  //     rhs); convergence_norm = rhs.NormL2();
  //     // std::cout << "\nConvergence Norm : " << convergence_norm <<
  //     std::endl;
  //   }
  //   if (convergence_norm > norm_goal && distance > tolerance) {
  //     std::cout << "DEFEAT" << std::endl;
  //   }
  // }
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
