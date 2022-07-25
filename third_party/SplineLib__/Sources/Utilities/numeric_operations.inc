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

template<typename Type>
constexpr Type GetEpsilon() {
  return Type{kEpsilonFactor<Type> * std::numeric_limits<Type>::epsilon()};
}

template<typename Type>
constexpr int GetPrecision() {
  return std::numeric_limits<Type>::digits10;
}

template<typename Type>
constexpr bool IsEqual(Type const &lhs, Type const &rhs, Type const &tolerance) {
#ifndef NDEBUG
  try {
    ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const &exception) {
    Throw(exception, "splinelib::sources::utilities::numeric_operations::IsEqual");
  }
#endif
  return std::abs(lhs - rhs) <= tolerance;
}

template<typename Type>
constexpr bool IsLess(Type const &lhs, Type const &rhs, Type const &tolerance) {
#ifndef NDEBUG
  try {
    ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const &exception) {
    Throw(exception, "splinelib::sources::utilities::numeric_operations::IsLess");
  }
#endif
  return IsEqual(lhs, rhs, tolerance) ? false : (lhs < rhs);
}

#ifndef NDEBUG
template<typename Type>
constexpr void ThrowIfToleranceIsNegative(Type const &tolerance) {
  tolerance >= Type{} ? void() : throw InvalidArgument("The provided tolerance " + std::to_string(tolerance) + " must "
                                                       "not be negative.");
}
#endif
