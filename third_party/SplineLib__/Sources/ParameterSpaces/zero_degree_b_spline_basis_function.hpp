/* Copyright (c) 2018–2021 SplineLib

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit
persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. */

#ifndef SOURCES_PARAMETERSPACES_ZERO_DEGREE_B_SPLINE_BASIS_FUNCTION_HPP_
#define SOURCES_PARAMETERSPACES_ZERO_DEGREE_B_SPLINE_BASIS_FUNCTION_HPP_

#include "Sources/ParameterSpaces/b_spline_basis_function.hpp"
#include "Sources/ParameterSpaces/knot_vector.hpp"
#include "Sources/Utilities/named_type.hpp"

namespace splinelib::sources::parameter_spaces {

// ZeroDegreeBSplineBasisFunctions N_{i,0} are indicator functions for the ith knot span.
//
// Example (see NURBS book Exa. 2.1):
//   KnotVector::Knot_ const k0_0{0.0}, k1_0{1.0};
//   KnotVector const kKnotVector{k0_0, k0_0, k0_0, k1_0, k1_0, k1_0};
//   ParametricCoordinate const k0_25{0.25};
//   ZeroDegreeBSplineBasisFunction const basis_function(kKnotVector, KnotSpan{2});
//   ZeroDegreeBSplineBasisFunction::Type_ const &evaluated = basis_function(k0_25);  // 1.0.
//                                         const &first_derivative = basis_function(k0_25, Derivative{1});  // 0.0.
class ZeroDegreeBSplineBasisFunction : public virtual BSplineBasisFunction {
 public:
  using Base_ = BSplineBasisFunction;
  using Type_ = Base_::Type_;

  ZeroDegreeBSplineBasisFunction() = default;
  ZeroDegreeBSplineBasisFunction(KnotVector const &knot_vector, KnotSpan const &start_of_support,
                                 Tolerance const &tolerance = kEpsilon);
  ZeroDegreeBSplineBasisFunction(ZeroDegreeBSplineBasisFunction const &other) = default;
  ZeroDegreeBSplineBasisFunction(ZeroDegreeBSplineBasisFunction &&other) noexcept = default;
  ZeroDegreeBSplineBasisFunction & operator=(ZeroDegreeBSplineBasisFunction const &rhs) = default;
  ZeroDegreeBSplineBasisFunction & operator=(ZeroDegreeBSplineBasisFunction &&rhs) noexcept = default;
  ~ZeroDegreeBSplineBasisFunction() override = default;

  // Comparison based on tolerance.
  friend bool IsEqual(ZeroDegreeBSplineBasisFunction const &lhs, ZeroDegreeBSplineBasisFunction const &rhs,
                      Tolerance const &tolerance);
  // Comparison based on numeric_operations::GetEpsilon<Tolerance>().
  friend bool operator==(ZeroDegreeBSplineBasisFunction const &lhs, ZeroDegreeBSplineBasisFunction const &rhs);

  Type_ operator()(ParametricCoordinate const &parametric_coordinate, Tolerance const &tolerance = kEpsilon) const
      override;
  Type_ operator()(ParametricCoordinate const &parametric_coordinate, Derivative const &derivative,
                   Tolerance const &tolerance = kEpsilon) const override;
};

bool IsEqual(ZeroDegreeBSplineBasisFunction const &lhs, ZeroDegreeBSplineBasisFunction const &rhs,
             Tolerance const &tolerance = kEpsilon);
bool operator==(ZeroDegreeBSplineBasisFunction const &lhs, ZeroDegreeBSplineBasisFunction const &rhs);

}  // namespace splinelib::sources::parameter_spaces

#endif  // SOURCES_PARAMETERSPACES_ZERO_DEGREE_B_SPLINE_BASIS_FUNCTION_HPP_
