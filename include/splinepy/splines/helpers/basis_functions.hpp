#pragma once

#include <numeric>

#include <bezman/src/utils/algorithms/recursive_combine.hpp>
#include <bezman/src/utils/fastbinomialcoefficient.hpp>

namespace splinepy::splines::helpers {

using bezman::utils::algorithms::RecursiveCombine;

/// Bezier basis - scalar io here
template<typename SplineType, typename QueryType, typename BasisValueType>
constexpr inline void BezierBasis(const SplineType& spline,
                                  const QueryType* para_coord,
                                  BasisValueType* basis) {
  // check if this is bezier spline. excuse me for ambiguity
  static_assert(!SplineType::kHasKnotVectors,
                "BezierBasis is only for bezier spline families.");

  // form query
  typename SplineType::ParametricCoordinate_ bz_query;
  std::copy_n(para_coord, SplineType::kParaDim, bz_query.begin());

  // query and copy output
  // no difference in call between rational and non-rational
  const auto bz_basis = spline.BasisFunctions(bz_query);
  std::copy(bz_basis.begin(), bz_basis.end(), basis);
}

/// Bezier basis derivative - scalar io here
template<typename SplineType,
         typename QueryType,
         typename OrderType,
         typename BasisValueType>
constexpr inline void BezierBasisDerivative(const SplineType& spline,
                                            const QueryType* para_coord,
                                            const OrderType* order,
                                            BasisValueType* basis_der) {
  // check if this is bezier spline. excuse me for ambiguity
  static_assert(!SplineType::kHasKnotVectors,
                "BezierBasis is only for bezier spline families.");

  // form query
  typename SplineType::ParametricCoordinate_ bz_query;
  typename SplineType::Derivative_ bz_order;
  for (int i{}; i < SplineType::kParaDim; ++i) {
    bz_query[i] =
        static_cast<typename SplineType::ParametricCoordinate_::value_type>(
            para_coord[i]);
    bz_order[i] =
        static_cast<typename SplineType::Derivative_::value_type>(order[i]);
  }

  // query and copy output
  // no difference in call between rational and non-rational
  const auto bz_basis_der =
      spline.BasisFunctionsDerivatives(bz_query, bz_order);
  std::copy(bz_basis_der.begin(), bz_basis_der.end(), basis_der);
}

// Bezier support - this is [0, n_support)
template<typename SplineType, typename QueryType, typename SupportType>
constexpr inline void BezierSupport(const SplineType& spline,
                                    const QueryType* para_coord,
                                    SupportType* support) {
  // check if this is bezier spline. excuse me for ambiguity
  static_assert(!SplineType::kHasKnotVectors,
                "BezierSupport is only for bezier spline families.");

  // fill
  std::iota(support, &support[spline.SplinepyNumberOfSupports()], 0);
}

/// BSpline Basis functions per dimension
template<typename SplineType, typename QueryType>
constexpr inline std::array<std::vector<double>, SplineType::kParaDim>
BSplineBasisPerParametricDimension(const SplineType& spline,
                                   const QueryType* para_coord) {

  static_assert(SplineType::kHasKnotVectors,
                "BSplineBasis is only for bspline families.");

  // prepare unique evals
  const auto& parameter_space = spline.GetParameterSpace();

  // prepare query
  typename SplineType::ParametricCoordinate_ sl_query;
  for (int i{}; i < SplineType::kParaDim; ++i) {
    sl_query[i] =
        typename SplineType::ScalarParametricCoordinate_{para_coord[i]};
  }

  return parameter_space.EvaluateBasisValuesPerDimension(sl_query);
}

/// BSpline Basis support
template<typename SplineType, typename QueryType>
inline std::vector<int> BSplineSupport(const SplineType& spline,
                                       const QueryType* para_coord) {
  static_assert(SplineType::kHasKnotVectors,
                "BSplineBasis is only for bspline families.");
  // prepare query
  typename SplineType::ParametricCoordinate_ sl_query;
  for (int i{}; i < SplineType::kParaDim; ++i) {
    sl_query[i] =
        typename SplineType::ScalarParametricCoordinate_{para_coord[i]};
  }

  // prepare properties for the loop
  const auto& parameter_space = spline.GetParameterSpace();
  const auto n_basis_per_dim =
      parameter_space.GetNumberOfNonZeroBasisFunctions();
  const auto support_start =
      parameter_space.FindFirstNonZeroBasisFunction(sl_query);
  auto support_offset = parameter_space.First();
  const int n_total_basis =
      std::accumulate(n_basis_per_dim.begin(),
                      n_basis_per_dim.end(),
                      typename decltype(n_basis_per_dim)::value_type{1},
                      std::multiplies{});

  // prepare return
  std::vector<int> supports;
  supports.reserve(n_total_basis);

  // fill supports.
  for (int i{}; i < n_total_basis; ++i) {
    supports.push_back(
        (support_start + support_offset.GetIndex()).GetIndex1d().Get());
    ++support_offset;
  }
  return supports;
}

template<typename SplineType, typename QueryType, typename SupportType>
inline void BSplineSupport(const SplineType& spline,
                           const QueryType* para_coord,
                           SupportType* out_support) {
  const auto support = BSplineSupport(spline, para_coord);
  std::copy(support.begin(), support.end(), out_support);
}

/// pure bspline basis. We recommend using BSplineBasis() instead of this
template<typename SplineType, typename QueryType>
constexpr inline std::vector<QueryType>
NonRationalBSplineBasis(const SplineType& spline, const QueryType* para_coord) {
  return RecursiveCombine(
      BSplineBasisPerParametricDimension(spline, para_coord));
}

/// nurbs basis. We recommend using BSplineBasis() instead of this
template<typename SplineType, typename QueryType, typename SupportType>
constexpr inline std::vector<QueryType>
RationalBSplineBasis(const SplineType& spline,
                     const QueryType* para_coord,
                     const SupportType& support) {
  static_assert(SplineType::kIsRational && SplineType::kHasKnotVectors,
                "RationalBSplineBasis is only applicable to NURBS.");
  const auto bspline_basis = NonRationalBSplineBasis(spline, para_coord);
  const auto& homogeneous_coords = spline.GetCoordinates();
  const auto n_basis = support.size();

  // prepare output
  std::vector<QueryType> rational_basis;
  rational_basis.reserve(n_basis);

  double W{0.};
  for (std::size_t i{}; i < n_basis; ++i) {
    // get weight
    const auto& w = homogeneous_coords[support[i]][SplineType::kDim];
    const auto N_times_w = bspline_basis[i] * w;

    W += N_times_w;
    rational_basis.push_back(N_times_w);
  }

  const auto W_inv = 1 / W;
  for (auto& rb : rational_basis) {
    rb *= W_inv;
  }

  return rational_basis;
}

/// nurbs basis. We recommend using BSplineBasis() instead of this
template<typename SplineType, typename QueryType>
constexpr inline std::vector<QueryType>
RationalBSplineBasis(const SplineType& spline, const QueryType* para_coord) {
  return RationalBSplineBasis(spline,
                              para_coord,
                              BSplineSupport(spline, para_coord));
}

/// BSpline Basis functions
template<typename SplineType, typename QueryType, typename BasisType>
constexpr inline void BSplineBasis(const SplineType& spline,
                                   const QueryType* para_coord,
                                   BasisType* basis) {

  if constexpr (SplineType::kIsRational) {
    const auto r_basis = RationalBSplineBasis(spline, para_coord);
    std::copy(r_basis.begin(), r_basis.end(), basis);
  } else {
    const auto nr_basis = NonRationalBSplineBasis(spline, para_coord);
    std::copy(nr_basis.begin(), nr_basis.end(), basis);
  }
}

/// BSpline Basis functions derivative per dimension
template<typename SplineType, typename QueryType, typename OrderType>
constexpr inline std::array<std::vector<double>, SplineType::kParaDim>
BSplineBasisDerivativePerParametricDimension(const SplineType& spline,
                                             const QueryType* para_coord,
                                             const OrderType* order) {

  static_assert(SplineType::kHasKnotVectors,
                "BSplineBasisDerivative is only for bspline families.");
  // prepare unique eval buckets
  const auto& parameter_space = spline.GetParameterSpace();

  // prepare query
  typename SplineType::ParametricCoordinate_ sl_query;
  typename SplineType::Derivative_ sl_order;
  for (int i{}; i < SplineType::kParaDim; ++i) {
    sl_query[i] =
        typename SplineType::ScalarParametricCoordinate_{para_coord[i]};
    sl_order[i] = typename SplineType::Derivative_::value_type{order[i]};
  }

  return parameter_space.EvaluateBasisDerivativeValuesPerDimension(sl_query,
                                                                   sl_order);
}

/// BSpline Basis functions der
template<typename SplineType, typename QueryType, typename OrderType>
constexpr inline std::vector<QueryType>
NonRationalBSplineBasisDerivative(const SplineType& spline,
                                  const QueryType* para_coord,
                                  const OrderType* order) {

  return RecursiveCombine(
      BSplineBasisDerivativePerParametricDimension(spline, para_coord, order));
}

/// adapted from bezman
template<typename SplineType, typename QueryType, typename OrderType>
constexpr inline std::vector<QueryType>
RationalBSplineBasisDerivative(const SplineType& spline,
                               const QueryType* para_coord,
                               const OrderType* order) {

  // we will do everything with OrderType
  constexpr auto para_dim = static_cast<OrderType>(SplineType::kParaDim);

  // prepare
  const auto& homogeneous_coords = spline.GetCoordinates();

  // weight getter shortcut
  auto get_weight_ = [&homogeneous_coords](const OrderType id) -> QueryType {
    return static_cast<QueryType>(homogeneous_coords[id][SplineType::kDim]);
  };

  // Define lambdas to switch between local indexing with
  // coordinate style indexing to global, scalar indexing
  // example: for nth (2,1,4) : (1,0,0)->1, (1,1,0)->4

  // Global (scalar) indexing to local index-system
  auto local_ids_ =
      [&](const OrderType req_id) -> std::array<OrderType, para_dim> {
    OrderType id{req_id};
    std::array<OrderType, para_dim> local_ids{};
    for (OrderType i_pd{}; i_pd < para_dim; ++i_pd) {
      if (order[i_pd] == 0)
        continue;
      local_ids[i_pd] = id % (order[i_pd] + 1);
      id -= local_ids[i_pd];
      id /= order[i_pd] + 1;
    }
    return local_ids;
  };

  // Local (coordinate-style) indexing to global
  auto global_ids_ = [&](const OrderType* req_derivs) -> OrderType {
    OrderType id{};
    OrderType offset{1};
    for (OrderType i_pd{}; i_pd < para_dim; ++i_pd) {
      // assert(req_derivs[i_pd] <= order[i_pd]);
      if (order[i_pd] == 0)
        continue;
      id += offset * (req_derivs[i_pd]);
      offset = order[i_pd] > 0 ? offset * (order[i_pd] + 1) : offset;
    }
    return id;
  };

  // Check if requested derivative is "subset" to current derivative
  auto is_not_subset_ =
      [&](const std::array<OrderType, para_dim>& req_derivs_max,
          const std::array<OrderType, para_dim>& req_derivs) -> bool {
    for (OrderType i_pd{}; i_pd < para_dim; ++i_pd) {
      if (req_derivs[i_pd] > req_derivs_max[i_pd])
        return true;
    }
    return false;
  };

  // Prepare supports
  const auto supports = BSplineSupport(spline, para_coord);
  // Initialize return type
  const OrderType number_of_derivs{global_ids_(order) + 1};
  const OrderType n_basis_functions{spline.SplinepyNumberOfSupports()};
  // Please remember that the first derivative is not used
  std::vector<std::vector<QueryType>> derivatives(number_of_derivs);
  std::vector<std::vector<QueryType>> A_derivatives(number_of_derivs);
  std::vector<QueryType> w_derivatives(number_of_derivs);

  // Fill all polynomial spline derivatives (and values for id=0)
  for (OrderType i_deriv{}; i_deriv < number_of_derivs; ++i_deriv) {
    const auto req_derivs = local_ids_(i_deriv);
    A_derivatives[i_deriv] =
        NonRationalBSplineBasisDerivative(spline,
                                          para_coord,
                                          req_derivs.data());
    for (OrderType i_basis{}; i_basis < n_basis_functions; ++i_basis) {
      const QueryType weight = get_weight_(supports[i_basis]);
      w_derivatives[i_deriv] += weight * A_derivatives[i_deriv][i_basis];
      A_derivatives[i_deriv][i_basis] *= weight;
    }
  }

  // Precompute inverse of weighted function
  const QueryType inv_w_fact = static_cast<QueryType>(1.) / w_derivatives[0];

  // Loop over all lower-order derivatives and assign derivatives-vector
  // Notation follows "The NURBS book" eq. 4.20 (extended for n-d splines)
  for (OrderType i_deriv{0}; i_deriv < number_of_derivs; ++i_deriv) {
    // Retrieve index-wise order of the derivative for current ID
    const auto derivative_order_indexwise_LHS = local_ids_(i_deriv);
    // Assign derivative of Numerator-function
    derivatives[i_deriv] = A_derivatives[i_deriv];
    // Substract all weighted lower-order functions
    for (OrderType j_deriv{1}; j_deriv <= i_deriv; ++j_deriv) {
      // Retrieve order of current index
      const auto derivative_order_indexwise_RHS = local_ids_(j_deriv);
      // Check only subsets  // Define lambdas to switch between local indexing
      // with coordinate style
      if (is_not_subset_(derivative_order_indexwise_LHS,
                         derivative_order_indexwise_RHS))
        continue;
      // Precompute Product of binomial coefficients
      // for now, use bezman's implementation - TODO: use SplineLib
      std::size_t binom_fact{1};
      for (OrderType i_pd{}; i_pd < para_dim; ++i_pd) {
        binom_fact *= bezman::utils::FastBinomialCoefficient::choose(
            derivative_order_indexwise_LHS[i_pd],
            derivative_order_indexwise_RHS[i_pd]);
      }
      // Substract low-order function
      for (OrderType i_basis{}; i_basis < n_basis_functions; ++i_basis) {
        derivatives[i_deriv][i_basis] -=
            static_cast<OrderType>(binom_fact) * w_derivatives[j_deriv]
            * derivatives[i_deriv - j_deriv][i_basis];
      }
    }
    // Finalize
    for (OrderType i_basis{}; i_basis < n_basis_functions; ++i_basis) {
      derivatives[i_deriv][i_basis] *= inv_w_fact;
    }
  }
  // Return last value
  return derivatives[number_of_derivs - 1];
}

template<typename SplineType,
         typename QueryType,
         typename OrderType,
         typename BasisType>
constexpr inline void BSplineBasisDerivative(const SplineType& spline,
                                             const QueryType* para_coord,
                                             const OrderType* order,
                                             BasisType* basis_der) {
  if constexpr (SplineType::kIsRational) {
    const auto r_basis =
        RationalBSplineBasisDerivative(spline, para_coord, order);
    std::copy(r_basis.begin(), r_basis.end(), basis_der);
  } else {
    const auto nr_basis =
        NonRationalBSplineBasisDerivative(spline, para_coord, order);
    std::copy(nr_basis.begin(), nr_basis.end(), basis_der);
  }
}

/// Nurbs Basis Functions der

} // namespace splinepy::splines::helpers
