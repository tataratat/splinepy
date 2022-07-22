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

#ifndef SOURCES_UTILITIES_NAMED_TYPE_HPP_
#define SOURCES_UTILITIES_NAMED_TYPE_HPP_

#include <cmath>
#include <functional>
#include <typeinfo>
#include <type_traits>
#include <utility>

#include "Sources/Utilities/error_handling.hpp"
#include "Sources/Utilities/numeric_operations.hpp"

namespace splinelib::sources::utilities {

template<typename Name, typename Type> class NamedType;

template<typename Name, typename Type>
constexpr NamedType<Name, Type> operator+(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs);
template<typename Name, typename Type>
constexpr NamedType<Name, Type> operator-(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs);
template<typename Name, typename Type>
constexpr NamedType<Name, Type> operator*(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs);
template<typename Name, typename Type>
constexpr NamedType<Name, Type> operator*(Type const &lhs, NamedType<Name, Type> const &rhs);
template<typename Name, typename Type>
constexpr NamedType<Name, Type> operator*(NamedType<Name, Type> const &lhs, Type const &rhs);
template<typename Name, typename Type>
constexpr NamedType<Name, Type> operator/(NamedType<Name, Type> const &dividend, NamedType<Name, Type> const &divisor);
template<typename Name, typename Type>
constexpr NamedType<Name, Type> operator/(Type const &dividend, NamedType<Name, Type> const &divisor);
template<typename Name, typename Type>
constexpr NamedType<Name, Type> operator/(NamedType<Name, Type> const &dividend, Type const &divisor);
template<typename Name, typename Type>
constexpr bool IsEqual(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs,
                       Type const &tolerance = numeric_operations::GetEpsilon<Type>());
template<typename Name, typename Type>
constexpr bool IsLessOrEqual(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs,
                             Type const &tolerance = numeric_operations::GetEpsilon<Type>());
template<typename Name, typename Type>
constexpr bool IsGreaterOrEqual(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs,
                                Type const &tolerance = numeric_operations::GetEpsilon<Type>());
template<typename Name, typename Type>
constexpr bool IsLess(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs,
                      Type const &tolerance = numeric_operations::GetEpsilon<Type>());
template<typename Name, typename Type>
constexpr bool IsGreater(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs,
                         Type const &tolerance = numeric_operations::GetEpsilon<Type>());
template<typename Name, typename Type>
constexpr bool operator==(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs);
template<typename Name, typename Type>
constexpr bool operator!=(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs);
template<typename Name, typename Type>
constexpr bool operator<=(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs);
template<typename Name, typename Type>
constexpr bool operator>=(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs);
template<typename Name, typename Type>
constexpr bool operator<(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs);
template<typename Name, typename Type>
constexpr bool operator>(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs);

// NamedTypes allow to check for semantics of variables (especially arguments) at compile time.  Actual connections
// between names and intrinsic types are provided by the definitions of NamedTypes placed at the end of this file.  In
// addition, they provide more reliable comparisons of floating point numbers.
//
// Example:
//   using Index = NamedType<struct IndexName, int>;
//   Index index{1};
//   Index::Type_ const value{(2 * ++index).Get()};  // Stores 4 in value.
//   using ParametricCoordinate = NamedType<struct ParametricCoordinateName, double>;
//   ParametricCoordinate const &parametric_coordinate = (ParametricCoordinate{1.0} + ParametricCoordinate{2.0});
//   parametric_coordinate -= ParametricCoordinate{3.0};  // Assigns ParametricCoordinate{} to parametric_coordinate.
//   ParametricCoordinate::Type_ const &epsilon = numeric_operations::GetEpsilon<ParametricCoordinate::Type_>();
//   parametric_coordinate == ParametricCoordinate{0.9 * epsilon};  // Equal with respect to tolerance = epsilon.
//   parametric_coordinate == ParametricCoordinate{1.1 * epsilon});  // Not equal with respect to tolerance = epsilon.
//   IsEqual(parametric_coordinate, ParametricCoordinate{1.1 * epsilon}, 1.2 * epsilon);  // Approximately equal.
//   constexpr bool const &kFalse = is_named_type<int>;
template<typename Name, typename Type>
class NamedType {
 public:
  using Name_ = Name;
  using Type_ = Type;

  constexpr NamedType();
  constexpr explicit NamedType(Type_ value);
  constexpr NamedType(NamedType const &other) = default;
  constexpr NamedType(NamedType &&other) noexcept = default;
  constexpr NamedType & operator=(NamedType const &rhs) = default;
  constexpr NamedType & operator=(NamedType &&rhs) noexcept = default;
  ~NamedType() = default;

  constexpr NamedType & operator+=(NamedType const &rhs);
  constexpr NamedType & operator-=(NamedType const &rhs);
  constexpr NamedType & operator++();
  constexpr NamedType & operator--();
  friend constexpr NamedType operator+<Name, Type>(NamedType const &lhs, NamedType const &rhs);
  friend constexpr NamedType operator-<Name, Type>(NamedType const &lhs, NamedType const &rhs);
  friend constexpr NamedType operator*<Name, Type>(NamedType const &lhs, NamedType const &rhs);
  friend constexpr NamedType operator*<Name, Type>(Type_ const &lhs, NamedType const &rhs);
  friend constexpr NamedType operator*<Name, Type>(NamedType const &lhs, Type_ const &rhs);
  friend constexpr NamedType operator/<Name, Type>(NamedType const &dividend, NamedType const &divisor);
  friend constexpr NamedType operator/<Name, Type>(Type_ const &dividend, NamedType const &divisor);
  friend constexpr NamedType operator/<Name, Type>(NamedType const &dividend, Type_ const &divisor);
  // Comparison based on given tolerance.
  friend constexpr bool IsEqual<Name, Type>(NamedType const &lhs, NamedType const &rhs, Type_ const &tolerance);
  friend constexpr bool IsLessOrEqual<Name, Type>(NamedType const &lhs, NamedType const &rhs, Type_ const &tolerance);
  friend constexpr bool IsGreaterOrEqual<Name, Type>(NamedType const &lhs, NamedType const &rhs,
                                                     Type_ const &tolerance);
  friend constexpr bool IsLess<Name, Type>(NamedType const &lhs, NamedType const &rhs, Type_ const &tolerance);
  friend constexpr bool IsGreater<Name, Type>(NamedType const &lhs, NamedType const &rhs, Type_ const &tolerance);
  // Comparison based on numeric_operations::GetEpsilon<Type>().
  friend constexpr bool operator==<Name, Type>(NamedType const &lhs, NamedType const &rhs);
  friend constexpr bool operator!=<Name, Type>(NamedType const &lhs, NamedType const &rhs);
  friend constexpr bool operator<=<Name, Type>(NamedType const &lhs, NamedType const &rhs);
  friend constexpr bool operator>=<Name, Type>(NamedType const &lhs, NamedType const &rhs);
  friend constexpr bool operator< <Name, Type>(NamedType const &lhs, NamedType const &rhs);
  friend constexpr bool operator> <Name, Type>(NamedType const &lhs, NamedType const &rhs);

  constexpr operator Type_() const {return value_;}
  constexpr operator Type_&() {return value_;}

  constexpr Type_ const & Get() const;
  constexpr Type_ & Get();

  static void ForEach(int const &first, int const &size,
                      std::function<void(NamedType<Name, int> const &)> const &function);

#ifndef NDEBUG
  static void ThrowIfNamedIntegerIsOutOfBounds(NamedType<Name, int> const &named_integer, int const &maximum_value);
#endif

 private:
  Type_ value_;
};

#include "Sources/Utilities/named_type.inc"

template<typename Type>
constexpr bool const is_named_type{false};
template<typename Name, typename Type>
constexpr bool const is_named_type<NamedType<Name, Type>>{true};

}  // namespace splinelib::sources::utilities

namespace splinelib {

// utilities
using Dimension = sources::utilities::NamedType<struct DimensionName, int>;
using Index = sources::utilities::NamedType<struct IndexName, int>;
using Length = sources::utilities::NamedType<struct LengthName, int>;
using Precision = sources::utilities::NamedType<struct PrecisionName, int>;

// parameter spaces
using Degree = sources::utilities::NamedType<struct DegreeName, int>;
using Derivative = sources::utilities::NamedType<struct DerivativeName, int>;
using KnotSpan = sources::utilities::NamedType<struct KnotSpanName, int>;
using Multiplicity = sources::utilities::NamedType<struct MultiplicityName, int>;
using ParametricCoordinate = sources::utilities::NamedType<struct ParametricCoordinateName, double>;

namespace sources::parameter_spaces {

using Type = ParametricCoordinate::Type_;
using BinomialRatio = Type;
using KnotRatio = Type;
using Tolerance = Type;

constexpr Multiplicity const kMultiplicity{1};
constexpr Precision const kPrecision{utilities::numeric_operations::GetPrecision<Type>()};
constexpr Tolerance const kEpsilon{utilities::numeric_operations::GetEpsilon<Tolerance>()};

}  // namespace sources::parameter_spaces

// vector spaces
using Coordinate = sources::utilities::NamedType<struct CoordinateName, sources::parameter_spaces::Type>;
using Weight = sources::utilities::NamedType<struct WeightName, Coordinate::Type_>;

namespace sources::vector_spaces {

using Type = Coordinate::Type_;
using Tolerance = Type;

constexpr Precision const kPrecision{utilities::numeric_operations::GetPrecision<Type>()};
constexpr Tolerance const kEpsilon{utilities::numeric_operations::GetEpsilon<Tolerance>()};

}  // namespace sources::vector_spaces

// splines
namespace sources::splines {

using Type = vector_spaces::Type;
using Tolerance = Type;

constexpr Multiplicity const kMultiplicity{parameter_spaces::kMultiplicity};
constexpr Precision const kPrecision{utilities::numeric_operations::GetPrecision<Type>()};
constexpr Tolerance const kEpsilon{utilities::numeric_operations::GetEpsilon<Tolerance>()};

}  // namespace sources::splines

// input output
namespace sources::input_output {

using Type = splines::Type;
using Tolerance = Type;

constexpr Precision const kPrecision{utilities::numeric_operations::GetPrecision<Type>()};
constexpr Tolerance const kEpsilon{utilities::numeric_operations::GetEpsilon<Tolerance>()};

}  // namespace sources::input_output

// models
namespace sources::models {

using Type = splines::Type;
using Tolerance = Type;

constexpr Precision const kPrecision{utilities::numeric_operations::GetPrecision<Type>()};
constexpr Tolerance const kEpsilon{utilities::numeric_operations::GetEpsilon<Tolerance>()};

}  // namespace sources::models

}  // namespace splinelib

#endif  // SOURCES_UTILITIES_NAMED_TYPE_HPP_
