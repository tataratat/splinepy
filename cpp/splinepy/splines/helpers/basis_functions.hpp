#pragma once

#include <numeric>

#include <bezman/src/utils/algorithms/recursive_combine.hpp>

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
  typename SplineType::ParameterSpace_::UniqueEvaluations_ uniq_evals;
  const auto& parameter_space = spline.GetParameterSpace();
  parameter_space.template InitializeUniqueEvaluations<false>(uniq_evals);

  // prepare query
  typename SplineType::ParametricCoordinate_ sl_query;
  for (int i{}; i < SplineType::kParaDim; ++i) {
    sl_query[i] =
        typename SplineType::ScalarParametricCoordinate_{para_coord[i]};
  }

  const auto n_basis_per_dim =
      parameter_space.GetNumberOfNonZeroBasisFunctions();
  const auto first_non_zeros =
      parameter_space.FindFirstNonZeroBasisFunction(sl_query);
  const auto& basis_functions = parameter_space.GetBasisFunctions();

  // prepare return
  std::array<std::vector<double>, SplineType::kParaDim> output_basis;

  // para dim loop
  for (int i{}; i < SplineType::kParaDim; ++i) {
    // things needed to compute
    const auto& non_zero_start = static_cast<int&>(first_non_zeros[i]);
    const auto& basis_per_dim = basis_functions[i];
    auto& uniq_evals_per_dim = uniq_evals[i];
    auto& sl_query_per_dim = sl_query[i];
    const auto n_basis = static_cast<int>(n_basis_per_dim[i]);
    // prepare output
    auto& basis_to_fill = output_basis[i];
    basis_to_fill.reserve(n_basis);
    // basis loop
    for (int j{}; j < static_cast<int>(n_basis); ++j) {
      // fill basis
      basis_to_fill.push_back(
          basis_per_dim[non_zero_start + j]->operator()(sl_query_per_dim,
                                                        uniq_evals_per_dim,
                                                        j));
    }
  }

  return output_basis;
}

/// BSpline Basis support
template<typename SplineType, typename QueryType>
inline std::vector<QueryType> BSplineSupport(const SplineType& spline,
                                             const QueryType* para_coord) {
  static_assert(SplineType::kHasKnotVectors,
                "BSplineBasis is only for bspline families.");
  // prepare query
  typename SplineType::ParametricCoordinate_ sl_query;
  for (int i{}; i < SplineType::kParaDim; ++i) {
    sl_query[i] =
        typename SplineType::ScalarParametricCoordinate_{para_coord[i]};
  }

  const auto& parameter_space = spline.GetParameterSpace();
  const auto n_basis_per_dim =
      parameter_space.GetNumberOfNonZeroBasisFunctions();
  auto support_index = parameter_space.FindFirstNonZeroBasisFunction(sl_query);
  const int n_total_basis =
      std::accumulate(n_basis_per_dim.begin(),
                      n_basis_per_dim.end(),
                      typename decltype(n_basis_per_dim)::value_type{1},
                      std::multiplies{})
          .Get();

  // prepare return
  std::vector<QueryType> supports;
  supports.reserve(n_total_basis);

  for (int i{}; i < n_total_basis; ++i) {
    supports.push_back(support_index.GetIndex1d().Get());
    ++support_index;
  }
  return supports;
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
  constexpr auto bspline_basis = NonRationalBSplineBasis(spline, para_coord);
  const auto& phys_space = spline.GetWeightedVectorSpace();
  const auto n_basis = support.size();

  // prepare output
  std::vector<QueryType> rational_basis;
  rational_basis.reserve(n_basis);

  double W{0.};
  for (std::size_t i{}; i < n_basis; ++i) {
    // get weight
    const auto& w = phys_space[support[i]][SplineType::kDim].Get();
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
template<typename SplineType, typename QueryType>
constexpr inline std::vector<QueryType>
BSplineBasis(const SplineType& spline, const QueryType* para_coord) {

  if constexpr (SplineType::kIsRational) {
    return RationalBSplineBasis(spline, para_coord);
  } else {
    return NonRationalBSplineBasis(spline, para_coord);
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
  typename SplineType::ParameterSpace_::UniqueDerivatives_ unique_derivatives;
  typename SplineType::ParameterSpace_::UniqueEvaluations_ unique_evaluations;
  typename SplineType::ParameterSpace_::IsTopLevelComputed_ top_levels_computed;
  const auto& parameter_space = spline.GetParameterSpace();

  // prepare query
  typename SplineType::ParametricCoordinate_ sl_query;
  typename SplineType::Derivative_ sl_order;
  for (int i{}; i < SplineType::kParaDim; ++i) {
    sl_query[i] =
        typename SplineType::ScalarParametricCoordinate_{para_coord[i]};
    sl_order[i] = typename SplineType::Derivative_::value_type{order[i]};
  }

  // initialize buckets
  parameter_space.template InitializeUniqueDerivativeContainers<false>(
      sl_order,
      top_levels_computed,
      unique_derivatives,
      unique_evaluations);

  // prepare loop data
  const auto n_basis_per_dim =
      parameter_space.GetNumberOfNonZeroBasisFunctions();
  const auto first_non_zeros =
      parameter_space.FindFirstNonZeroBasisFunction(sl_query);
  const auto& basis_functions = parameter_space.GetBasisFunctions();

  // prepare return
  std::array<std::vector<double>, SplineType::kParaDim> output_ders;

  // para dim loop
  for (int i{}; i < SplineType::kParaDim; ++i) {
    // things needed to compute
    const auto& non_zero_start = static_cast<int&>(first_non_zeros[i]);
    const auto& basis_function = basis_functions[i];
    auto& unique_evaluation = unique_evaluations[i];
    auto& unique_derivative = unique_derivatives[i];
    auto& top_level_computed = top_levels_computed[i];
    const auto& sl_query_entry = sl_query[i];
    const auto& sl_order_entry = sl_order[i];
    const auto n_basis = static_cast<int>(n_basis_per_dim[i]);
    // prepare output
    auto& basis_to_fill = output_ders[i];
    basis_to_fill.reserve(n_basis);

    if (order[i] != 0) {
      // der query
      for (int j{}; j < static_cast<int>(n_basis); ++j) {
        basis_to_fill.push_back(
            basis_function[non_zero_start + j]->operator()(sl_query_entry,
                                                           sl_order_entry,
                                                           unique_derivative,
                                                           unique_evaluation,
                                                           top_level_computed,
                                                           j));
      }
    } else {
      // this is just eval
      for (int j{}; j < static_cast<int>(n_basis); ++j) {
        // fill basis
        basis_to_fill.push_back(
            basis_function[non_zero_start + j]->operator()(sl_query_entry,
                                                           unique_evaluation,
                                                           j));
      }
    }
  }

  return output_ders;
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

/// Nurbs Basis Functions der

} // namespace splinepy::splines::helpers
