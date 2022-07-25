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

#include "Sources/ParameterSpaces/zero_degree_b_spline_basis_function.hpp"

#include "Sources/ParameterSpaces/knot_vector.hpp"
#include "Sources/Utilities/error_handling.hpp"
#include "Sources/Utilities/numeric_operations.hpp"

namespace splinelib::sources::parameter_spaces {

ZeroDegreeBSplineBasisFunction::ZeroDegreeBSplineBasisFunction(KnotVector const &knot_vector,
    KnotSpan const &start_of_support, Tolerance const &tolerance) : BSplineBasisFunction(knot_vector, start_of_support,
                                                                                         Degree{}, tolerance) {}

bool IsEqual(ZeroDegreeBSplineBasisFunction const &lhs, ZeroDegreeBSplineBasisFunction const &rhs,
             Tolerance const &tolerance) {
  using Base = ZeroDegreeBSplineBasisFunction::Base_;

#ifndef NDEBUG
  try {
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const &exception) {
    Throw(exception, "splinelib::sources::parameter_spaces::IsEqual::ZeroDegreeBSplineBasisFunction");
  }
#endif
  return IsEqual(static_cast<Base const &>(lhs), static_cast<Base const &>(rhs), tolerance);
}

bool operator==(ZeroDegreeBSplineBasisFunction const &lhs, ZeroDegreeBSplineBasisFunction const &rhs) {
  return IsEqual(lhs, rhs);
}

ZeroDegreeBSplineBasisFunction::Type_
ZeroDegreeBSplineBasisFunction::operator()(ParametricCoordinate const &parametric_coordinate,
                                           Tolerance const &tolerance) const {
#ifndef NDEBUG
  try {
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const &exception) {
    Throw(exception, "splinelib::sources::parameter_spaces::ZeroDegreeBSplineBasisFunction::operator()");
  }
#endif
  return (IsInSupport(parametric_coordinate, tolerance) ? Type_{1.0} : Type_{});
}

ZeroDegreeBSplineBasisFunction::Type_
ZeroDegreeBSplineBasisFunction::operator()(ParametricCoordinate const &parametric_coordinate,
                                           Derivative const &derivative, Tolerance const &tolerance) const {
#ifndef NDEBUG
  try {
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const &exception) {
    Throw(exception, "splinelib::sources::parameter_spaces::ZeroDegreeBSplineBasisFunction::operator()");
  }
#endif
  return (derivative == Derivative{} ? operator()(parametric_coordinate, tolerance) : Type_{});
}

}  // namespace splinelib::sources::parameter_spaces
