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
template<typename SplineType, typename QueryType, typename IntType>
void ScalarTypeElevateDegree(SplineType& spline,
                             const QueryType query,
                             const IntType multiplicity = 1) {
  using Dim = typename SplineType::Dimension_;

  if constexpr (SplineType::kHasKnotVectors) {
    // BSplineLib can take multiplicity
    spline.ElevateDegree(static_cast<Dim>(query), multiplicity);
  } else {
    for (IntType i{}; i < multiplicity; ++i) {
      spline.ElevateDegree(static_cast<Dim>(query));
    }
  }
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

/// @brief Inserts knot at given location with given multiplicity. It returns
/// number of successful multiplicity
///
/// @tparam SplineType
/// @tparam QueryDimType
/// @tparam QueryType
/// @tparam IntType
/// @param spline
/// @param query_dim
/// @param query
/// @param multiplicity
/// @return
template<typename SplineType,
         typename QueryDimType,
         typename QueryType,
         typename IntType>
int ScalarTypeInsertKnot(SplineType& spline,
                         const QueryDimType query_dim,
                         const QueryType query,
                         const IntType multiplicity) {

  using Dim = typename SplineType::Dimension_;
  using Knot = typename SplineType::Knot_;

  // since BSpline::InsertKnot is void func, we count knots before and after to
  // see if it was successful
  const auto& knot_vector =
      *spline.GetParameterSpace().GetKnotVectors()[query_dim];

  // Check if request is valid before BSplineLib
  // see discussions in PR #297
  // BSplineLib now also has runtime bound check, which raises error.
  const auto para_bounds = GetParametricBounds(spline);
  if ((para_bounds[0][query_dim] >= query)
      || (para_bounds[1][query_dim] <= query)) {
    return 0;
  }

  const auto n_knots_before = knot_vector.GetSize();

  spline.InsertKnot(Dim{query_dim}, Knot{query}, multiplicity);

  const auto n_knots_after = knot_vector.GetSize();

  return (n_knots_after - n_knots_before);
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
