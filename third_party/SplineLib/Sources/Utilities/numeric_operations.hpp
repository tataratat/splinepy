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

#ifndef SOURCES_UTILITIES_NUMERIC_OPERATIONS_HPP_
#define SOURCES_UTILITIES_NUMERIC_OPERATIONS_HPP_

#include <cmath>
#include <limits>

#include "Sources/Utilities/error_handling.hpp"

// Numeric operations such as 1.) getting the default tolerance and default precision for a Type and 2.) ordering Type
// values with respect to some tolerance.  Some Type's default tolerance is based on its machine epsilon (cf.
// GetEpsilon<Type>).  Note that the machine epsilon is only meaningful if std::numeric_limits<Type>::is_integer is
// false.
//
// Example:
//   double const &default_tolerance = GetEpsilon<double>();  // Approximately 2.22045e-14 (depending on machine).
//   int const &default_precision = GetPrecision<double>();  // The default precision is 16 digits.
//   IsEqual(0.0, 2e-14);  // Returns true as |0.0 - 2e-14| is less than or equal to the default tolerance.
//   IsEqual(0.0, 3e-14);  // Returns false as |0.0 - 3e-14| is greater than the default tolerance.
//   IsEqual(0.0, 3e-14, 3e-14);  // Returns true as |0.0 - 3e-14| is less than or equal to the specified tolerance.
namespace splinelib::sources::utilities::numeric_operations {

template<typename Type>
constexpr Type const kEpsilonFactor{100};

template<typename Type>
constexpr Type GetEpsilon();
template<typename Type>
constexpr int GetPrecision();
template<typename Type>
constexpr bool IsEqual(Type const &lhs, Type const &rhs, Type const &tolerance = GetEpsilon<Type>());
template<typename Type>
constexpr bool IsLess(Type const &lhs, Type const &rhs, Type const &tolerance = GetEpsilon<Type>());

#ifndef NDEBUG
template<typename Type>
constexpr void ThrowIfToleranceIsNegative(Type const &tolerance);
#endif

#include "Sources/Utilities/numeric_operations.inc"

}  // namespace splinelib::sources::utilities::numeric_operations

#endif  // SOURCES_UTILITIES_NUMERIC_OPERATIONS_HPP_
