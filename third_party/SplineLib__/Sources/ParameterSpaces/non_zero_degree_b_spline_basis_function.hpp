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

#ifndef SOURCES_PARAMETERSPACES_NON_ZERO_DEGREE_B_SPLINE_BASIS_FUNCTION_HPP_
#define SOURCES_PARAMETERSPACES_NON_ZERO_DEGREE_B_SPLINE_BASIS_FUNCTION_HPP_

#include "Sources/ParameterSpaces/b_spline_basis_function.hpp"
#include "Sources/ParameterSpaces/knot_vector.hpp"
#include "Sources/Utilities/named_type.hpp"
#include "Sources/Utilities/std_container_operations.hpp"

namespace splinelib::sources::parameter_spaces {

// NonZeroDegreeBSplineBasisFunctions N_{i,p} are piecewise polynomial functions of degree p > 0.  They are linear
// combinations of the basis functions N_{i,p-1} and N_{i+1,p-1} (see NURBS book Eq. (2.5)).  This recursion formula is
// implemented as a data structure. Note that 1/0 := 0.
//
// Example (see NURBS book Exa. 2.1):
//   KnotVector::Knot_ const k0_0{0.0}, k1_0{1.0};
//   KnotVector const kKnotVector{k0_0, k0_0, k0_0, k1_0, k1_0, k1_0};
//   ParametricCoordinate const k0_25{0.25};
//   NonZeroDegreeBSplineBasisFunction const basis_function(kKnotVector, KnotSpan{}, Degree{2});
//   NonZeroDegreeBSplineBasisFunction::Type_ const &evaluated = basis_function(k0_25);  // 0.5625.
//                                            const &first_derivative = basis_function(k0_25, Derivative{1});  // -1.5.
class NonZeroDegreeBSplineBasisFunction : public virtual BSplineBasisFunction {
 public:
  using Base_ = BSplineBasisFunction;
  using Type_ = Base_::Type_;

  NonZeroDegreeBSplineBasisFunction() = default;
  NonZeroDegreeBSplineBasisFunction(KnotVector const &knot_vector, KnotSpan const &start_of_support, Degree degree,
                                    Tolerance const &tolerance = kEpsilon);
  NonZeroDegreeBSplineBasisFunction(NonZeroDegreeBSplineBasisFunction const &other) = default;
  NonZeroDegreeBSplineBasisFunction(NonZeroDegreeBSplineBasisFunction &&other) noexcept = default;
  NonZeroDegreeBSplineBasisFunction & operator=(NonZeroDegreeBSplineBasisFunction const &rhs) = default;
  NonZeroDegreeBSplineBasisFunction & operator=(NonZeroDegreeBSplineBasisFunction &&rhs) noexcept = default;
  ~NonZeroDegreeBSplineBasisFunction() override = default;

  // Comparison based on tolerance.
  friend bool IsEqual(NonZeroDegreeBSplineBasisFunction const &lhs, NonZeroDegreeBSplineBasisFunction const &rhs,
                      Tolerance const &tolerance);
  // Comparison based on numeric_operations::GetEpsilon<Tolerance>().
  friend bool operator==(NonZeroDegreeBSplineBasisFunction const &lhs, NonZeroDegreeBSplineBasisFunction const &rhs);
  Type_ operator()(ParametricCoordinate const &parametric_coordinate, Tolerance const &tolerance = kEpsilon) const
      override;
  Type_ operator()(ParametricCoordinate const &parametric_coordinate, Derivative const &derivative,
                   Tolerance const &tolerance = kEpsilon) const override;

 protected:
  Type_ left_denominator_inverse_, right_denominator_inverse_, left_quotient_derivative_, right_quotient_derivative_;
  // Shared instead of unique pointers as BSplineBasisFunctions cannot be modified after their construction.
  SharedPointer<Base_> left_lower_degree_basis_function_, right_lower_degree_basis_function_;

 private:
  Type_ InvertPotentialZero(ParametricCoordinate const &potential_zero, Tolerance const &tolerance) const;  // 1/0 := 0.
};

bool IsEqual(NonZeroDegreeBSplineBasisFunction const &lhs, NonZeroDegreeBSplineBasisFunction const &rhs,
             Tolerance const &tolerance = kEpsilon);
bool operator==(NonZeroDegreeBSplineBasisFunction const &lhs, NonZeroDegreeBSplineBasisFunction const &rhs);

}  // namespace splinelib::sources::parameter_spaces

#endif  // SOURCES_PARAMETERSPACES_NON_ZERO_DEGREE_B_SPLINE_BASIS_FUNCTION_HPP_
