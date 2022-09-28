#pragma once

#include <algorithm>

namespace splinepy::splines::helpers {
/// SplineLib spline evaluation (single query).
template<typename SplineType,
         typename QueryType,
         typename OutputType>
void RawPtrEvaluate(SplineType& spline,
                    QueryType* query,
                    OutputType* output) {
  using Query = typename SplineType::ParametricCoordinate_;
  using QueryValueType = typename Query::value_type;
  Query core_query;

  for (std::size_t i{}; i < SplineType::kParaDim; ++i) {
    core_query[i] = QueryValueType{query[i]};
  }

  // @jzwar maybe somesort of if constexpr or maybe a new one for Bezier
  const auto core_evaluated = spline(core_query);

  for (std::size_t i{}; i < SplineType::kDim; ++i) {
    output[i] = static_cast<OutputType>(core_evaluated[i]);
  }
}

} // namespace splinepy::splines::helpers
