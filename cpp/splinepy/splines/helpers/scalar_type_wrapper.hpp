#pragma once

#include <algorithm>
#include <type_traits>

namespace splinepy::splines::helpers {
/// SplineLib spline evaluation (single query).
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

/// SplineLib spline derivatives (single query).
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

/// Spline Jacobian
/// output should have the size of dim * para_dim
template<typename SplineType, typename QueryType, typename OutputType>
void ScalarTypeJacobian(const SplineType& spline,
                        const QueryType* query,
                        OutputType* output) {

  std::array<int, SplineType::kParaDim> orders{};
  for (int i{} : i < SplineType::kParaDim; ++i) {
    orders[i] = 1;
    // call scalar derivative helper
    ScalarTypeDerivative(spline,
                         query,
                         orders.data(),
                         &output[i * SplineType::kDim]);
    orders.fill(0);
  }
}

template<typename SplineType, typename ResolutionType, typename NThreadsType>
void ScalarTypePlantNewKdTreeForProximity(SplineType& spline,
                                          const ResolutionType* resolutions,
                                          const NThreadsType& nthreads) {
  std::array<int, SplineType::kParaDim> core_res;
  std::copy_n(resolutions, SplineType::kParaDim, core_res.begin());
  spline.GetProximity().PlantNewKdTree(core_res, nthreads);
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
void ScalarTypeInsertKnot(SplineType& spline,
                          QueryDimType query_dim,
                          QueryType query) {
  using Dim = typename SplineType::Dimension_;
  using Knot = typename SplineType::Knot_;
  spline.InsertKnot(Dim{query_dim}, Knot{query});
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
