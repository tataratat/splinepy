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

#include "Sources/ParameterSpaces/non_zero_degree_b_spline_basis_function.hpp"

#include <utility>

#include "Sources/ParameterSpaces/knot_vector.hpp"
#include "Sources/Utilities/error_handling.hpp"
#include "Sources/Utilities/numeric_operations.hpp"

namespace splinelib::sources::parameter_spaces {

NonZeroDegreeBSplineBasisFunction::NonZeroDegreeBSplineBasisFunction(KnotVector const &knot_vector,
    KnotSpan const &start_of_support, Degree degree, Tolerance const &tolerance) : Base_(knot_vector, start_of_support,
        std::move(degree), tolerance), left_lower_degree_basis_function_{CreateDynamic(knot_vector, start_of_support,
            degree_ - Degree{1}, tolerance)}, right_lower_degree_basis_function_{CreateDynamic(knot_vector,
                start_of_support + KnotSpan{1}, degree_ - Degree{1}, tolerance)} {
    Index const start_index{start_of_support.Get()};
    Degree::Type_ const &degree_value = degree_.Get();
    left_denominator_inverse_ = InvertPotentialZero(knot_vector[start_index + Index{degree_value}] - start_knot_,
                                                    tolerance);
    right_denominator_inverse_ = InvertPotentialZero(end_knot_ - knot_vector[start_index + Index{1}], tolerance);
    left_quotient_derivative_ = (degree_value * left_denominator_inverse_);
    right_quotient_derivative_ = (degree_value * right_denominator_inverse_);
}

bool IsEqual(NonZeroDegreeBSplineBasisFunction const &lhs, NonZeroDegreeBSplineBasisFunction const &rhs,
             Tolerance const &tolerance) {
  using Base = NonZeroDegreeBSplineBasisFunction::Base_;
  using utilities::numeric_operations::IsEqual;

#ifndef NDEBUG
  try {
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const &exception) {
    Throw(exception, "splinelib::sources::parameter_spaces::IsEqual::NonZeroDegreeBSplineBasisFunction");
  }
#endif
  return (IsEqual(lhs.left_denominator_inverse_, rhs.left_denominator_inverse_, tolerance) &&
          IsEqual(lhs.right_denominator_inverse_, rhs.right_denominator_inverse_, tolerance) &&
          IsEqual(lhs.left_quotient_derivative_, rhs.left_quotient_derivative_, tolerance) &&
          IsEqual(lhs.right_quotient_derivative_, rhs.right_quotient_derivative_, tolerance) &&
          IsEqual(*lhs.left_lower_degree_basis_function_, *rhs.left_lower_degree_basis_function_, tolerance) &&
          IsEqual(*lhs.right_lower_degree_basis_function_, *rhs.right_lower_degree_basis_function_, tolerance) &&
          IsEqual(static_cast<Base const &>(lhs), static_cast<Base const &>(rhs), tolerance));
}

bool operator==(NonZeroDegreeBSplineBasisFunction const &lhs, NonZeroDegreeBSplineBasisFunction const &rhs) {
  return IsEqual(lhs, rhs);
}

// Recurrence formula due to DeBoor, Cox, and Mansfield (see NURBS book Eq. (2.5)).
NonZeroDegreeBSplineBasisFunction::Type_
NonZeroDegreeBSplineBasisFunction::operator()(ParametricCoordinate const &parametric_coordinate,
                                              Tolerance const &tolerance) const {
#ifndef NDEBUG
  try {
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const &exception) {
    Throw(exception, "splinelib::sources::parameter_spaces::NonZeroDegreeBSplineBasisFunction::operator()");
  }
#endif
  return IsInSupport(parametric_coordinate, tolerance) ? (((parametric_coordinate - start_knot_).Get() *
      left_denominator_inverse_ * (*left_lower_degree_basis_function_)(parametric_coordinate, tolerance)) +
          ((end_knot_ - parametric_coordinate).Get() * right_denominator_inverse_ *
           (*right_lower_degree_basis_function_)(parametric_coordinate, tolerance))) : Type_{};
}

// Based on recurrence formula due to DeBoor, Cox, and Mansfield (see NURBS book Eq. (2.9)).
NonZeroDegreeBSplineBasisFunction::Type_
NonZeroDegreeBSplineBasisFunction::operator()(ParametricCoordinate const &parametric_coordinate,
                                              Derivative const &derivative, Tolerance const &tolerance) const {
#ifndef NDEBUG
  try {
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const &exception) {
    Throw(exception, "splinelib::sources::parameter_spaces::NonZeroDegreeBSplineBasisFunction::operator()");
  }
#endif
  if (derivative != Derivative{}) {
    if (IsInSupport(parametric_coordinate, tolerance)) {
      Derivative const lower_derivative = (derivative - Derivative{1});
      return ((left_quotient_derivative_ * (*left_lower_degree_basis_function_)(parametric_coordinate, lower_derivative,
                   tolerance)) - (right_quotient_derivative_ *
                       (*right_lower_degree_basis_function_)(parametric_coordinate, lower_derivative, tolerance)));
    } else {
      return Type_{};
    }
  } else {
    return operator()(parametric_coordinate, tolerance);
  }
}

NonZeroDegreeBSplineBasisFunction::Type_
NonZeroDegreeBSplineBasisFunction::InvertPotentialZero(ParametricCoordinate const &potential_zero,
                                                       Tolerance const &tolerance) const {
  return IsEqual(potential_zero, ParametricCoordinate{}, tolerance) ? Type{} : Type_{1.0 / potential_zero.Get()};
}

}  // namespace splinelib::sources::parameter_spaces
