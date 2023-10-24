#pragma once

#include <algorithm>
#include <type_traits>

#include <splinepy/utils/arrays.hpp>
#include <splinepy/utils/default_initialization_allocator.hpp>

namespace splinepy::splines::helpers {
/// bezman spline evaluation (single query).
template<typename SplineType, typename QueryType, typename OutputType>
void ScalarTypeEvaluate(const SplineType& spline,
                        const QueryType* query,
                        OutputType* output) {
  using Query = typename SplineType::ParametricCoordinate_;
  using QueryValueType = typename Query::value_type;
  Query core_query;

  // form query
  for (std::size_t i{}; i < SplineType::kParaDim; ++i) {
    core_query[i] = QueryValueType{query[i]};
  }

  // query
  const auto core_evaluated = spline(core_query);

  // fill output
  if constexpr (std::is_scalar<decltype(core_evaluated)>::value) {
    output[0] = static_cast<OutputType>(core_evaluated);
  } else {
    for (std::size_t i{}; i < SplineType::kDim; ++i) {
      output[i] = static_cast<OutputType>(core_evaluated[i]);
    }
  }
}

/// @brief Computes parametric coordinate of boundary centers
/// @tparam SplineType
/// @tparam OutputType
/// @param spline
/// @param output should have size of 2 * para_dim * para_dim
template<typename SplineType, typename OutputType>
void ScalarTypeBoundaryCenters(const SplineType& spline, OutputType* output) {
  using DoubleVector = splinepy::utils::DefaultInitializationVector<double>;

  // Prepare inputs
  const int para_dim = spline.SplinepyParaDim();

  // get parametric bounds
  // They are given back in the order [min_0, min_1,...,max_0, max_1...]
  DoubleVector bounds_vector(2 * para_dim);
  double* bounds = bounds_vector.data();
  spline.SplinepyParametricBounds(bounds);

  // Set parametric coordinate queries
  for (int i{}; i < para_dim; ++i) {
    for (int j{}; j < para_dim; ++j) {
      if (i == j) {
        output[2 * i * para_dim + j] = bounds[j];
        output[(2 * i + 1) * para_dim + j] = bounds[j + para_dim];
      } else {
        const auto q = .5 * (bounds[j] + bounds[j + para_dim]);
        output[2 * i * para_dim + j] = q;
        output[(2 * i + 1) * para_dim + j] = q;
      }
    }
  }
}

/// @brief Evaluate Splines at boundary face centers
/// output should have size of 2 * para_dim * dim
template<typename SplineType, typename OutputType>
void ScalarTypeEvaluateBoundaryCenters(const SplineType& spline,
                                       OutputType* output) {

  using DoubleVector = splinepy::utils::DefaultInitializationVector<double>;

  // Prepare inputs
  const int para_dim = spline.SplinepyParaDim();
  const int dim = spline.SplinepyDim();
  const int n_faces = 2 * para_dim;

  DoubleVector queries_vector(n_faces * para_dim);
  double* queries = queries_vector.data();

  // compute queries
  ScalarTypeBoundaryCenters(spline, queries);

  // Evaluate
  for (int i{}; i < n_faces; ++i) {
    spline.SplinepyEvaluate(&queries[i * para_dim], &output[i * dim]);
  }
}

/// bezman spline derivatives (single query).
template<typename SplineType,
         typename QueryType,
         typename OrderType,
         typename OutputType>
void ScalarTypeDerivative(const SplineType& spline,
                          const QueryType* query,
                          const OrderType* order,
                          OutputType* output) {
  using Query = typename SplineType::ParametricCoordinate_;
  using QueryValueType = typename Query::value_type;
  using Order = typename SplineType::Derivative_;
  using OrderValueType = typename Order::value_type;
  Query core_query;
  Order core_order;

  // form query
  for (std::size_t i{}; i < SplineType::kParaDim; ++i) {
    core_query[i] = QueryValueType{query[i]};
    core_order[i] = static_cast<OrderValueType>(order[i]);
  }

  // query
  const auto core_derived = spline(core_query, core_order);

  // fill output
  if constexpr (std::is_scalar<decltype(core_derived)>::value) {
    output[0] = static_cast<OutputType>(core_derived);
  } else {
    for (std::size_t i{}; i < SplineType::kDim; ++i) {
      output[i] = static_cast<OutputType>(core_derived[i]);
    }
  }
}

/**
 * @brief inserts jacobian into requested pointer for single query
 *
 * @tparam SplineType Type of core spline
 * @tparam QueryType Spline-Type dependent query type
 * @tparam OutputType mostly double*
 * @param spline spline for evaluation (core)
 * @param query position for evaluation
 * @param output
 */
template<typename SplineType, typename QueryType, typename OutputType>
void ScalarTypeJacobian(const SplineType& spline,
                        const QueryType* query,
                        OutputType* output) {
  using RealArray = splinepy::utils::Array<double>;
  using IntArray = splinepy::utils::Array<int>;
  using RealArray2D = splinepy::utils::Array<double, 2>;

  const int dim = spline.SplinepyDim();
  // create derivative query, temp result holder, view on output
  IntArray der_query(SplineType::kParaDim);
  der_query.Fill(0);
  RealArray der_result(dim);
  RealArray2D output_view(output, dim, SplineType::kParaDim);

  // para_dim loop
  for (int i{}; i < SplineType::kParaDim; ++i) {
    // prepare eye query
    ++der_query[i];

    // derivative query
    spline.SplinepyDerivative(query, der_query.data(), der_result.data());

    // transposed fill
    for (int j{}; j < dim; ++j) {
      output_view(j, i) = der_result[j];
    }

    // reset eye query to zero
    --der_query[i];
  }
}

/// single degree elevation.
template<typename SplineType, typename QueryType>
void ScalarTypeElevateDegree(SplineType& spline, const QueryType query) {
  using Dim = typename SplineType::Dimension_;
  spline.ElevateDegree(static_cast<Dim>(query));
}

/// single degree reduction
template<typename SplineType, typename QueryType, typename ToleranceType>
bool ScalarTypeReduceDegree(SplineType& spline,
                            const QueryType query,
                            const ToleranceType tolerance) {
  using Dim = typename SplineType::Dimension_;
  using Tol = typename SplineType::Tolerance_;
  return spline.ReduceDegree(Dim{query}, Tol{tolerance});
}

/// single knot insertion
template<typename SplineType, typename QueryDimType, typename QueryType>
bool ScalarTypeInsertKnot(SplineType& spline,
                          QueryDimType query_dim,
                          QueryType query) {

  using Dim = typename SplineType::Dimension_;
  using Knot = typename SplineType::Knot_;

  // since BSpline::InsertKnot is void func, we count knots before and after to
  // see if it was successful
  const auto& knot_vector =
      *spline.GetParameterSpace().GetKnotVectors()[query_dim];
  const auto n_knots_before = knot_vector.GetSize();

  spline.InsertKnot(Dim{query_dim}, Knot{query});

  const auto n_knots_after = knot_vector.GetSize();

  return (n_knots_after > n_knots_before);
}

/// single knot removal
template<typename SplineType,
         typename QueryDimType,
         typename QueryType,
         typename ToleranceType>
bool ScalarTypeRemoveKnot(SplineType& spline,
                          QueryDimType query_dim,
                          QueryType query,
                          ToleranceType tolerance) {
  using Dim = typename SplineType::Dimension_;
  using Knot = typename SplineType::Knot_;
  using Tol = typename SplineType::Tolerance_;

  const auto multiplicity =
      spline.RemoveKnot(Dim{query_dim}, Knot{query}, Tol{tolerance});

  // very confusing syntax, let's see if this is correct
  // TODO: is this correct?
  if (static_cast<int>(multiplicity) == 0) {
    return false;
  }
  return true;
}

} // namespace splinepy::splines::helpers
