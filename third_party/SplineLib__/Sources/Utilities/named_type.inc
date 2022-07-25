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

template<typename Name, typename Type>
constexpr NamedType<Name, Type>::NamedType() : value_{} {}

template<typename Name, typename Type>
constexpr NamedType<Name, Type>::NamedType(Type_ value) : value_{std::move(value)} {
#ifndef NDEBUG
  if constexpr (std::is_same_v<Type_, int>) {
    value_ >= Type_{} ? void() : throw DomainError("Negative named integer (" + std::to_string(value_) + ") "
                                                   "encountered, but named integers must not be negative.");
  }
#endif
}

template<typename Name, typename Type>
constexpr NamedType<Name, Type> & NamedType<Name, Type>::operator+=(NamedType const &rhs) {
  return (*this = (*this + rhs));
}

template<typename Name, typename Type>
constexpr NamedType<Name, Type> & NamedType<Name, Type>::operator-=(NamedType const &rhs) {
  return (*this = (*this - rhs));
}

template<typename Name, typename Type>
constexpr NamedType<Name, Type> & NamedType<Name, Type>::operator++() {
  ++value_;
  return *this;
}

template<typename Name, typename Type>
constexpr NamedType<Name, Type> & NamedType<Name, Type>::operator--() {
  --value_;
  return *this;
}

template<typename Name, typename Type>
constexpr NamedType<Name, Type> operator+(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs) {
  return NamedType<Name, Type>{lhs.value_ + rhs.value_};
}

template<typename Name, typename Type>
constexpr NamedType<Name, Type> operator-(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs) {
  return NamedType<Name, Type>{lhs.value_ - rhs.value_};
}

template<typename Name, typename Type>
constexpr NamedType<Name, Type> operator*(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs) {
  return NamedType<Name, Type>{lhs.value_ * rhs.value_};
}

template<typename Name, typename Type>
constexpr NamedType<Name, Type> operator*(Type const &lhs, NamedType<Name, Type> const &rhs) {
  return NamedType<Name, Type>{lhs * rhs.value_};
}

template<typename Name, typename Type>
constexpr NamedType<Name, Type> operator*(NamedType<Name, Type> const &lhs, Type const &rhs) {
  return NamedType<Name, Type>{lhs.value_ * rhs};
}

template<typename Name, typename Type>
constexpr NamedType<Name, Type> operator/(NamedType<Name, Type> const &dividend, NamedType<Name, Type> const &divisor) {
  return (dividend / divisor.value_);
}

template<typename Name, typename Type>
constexpr NamedType<Name, Type> operator/(Type const &dividend, NamedType<Name, Type> const &divisor) {
  using NamedType = NamedType<Name, Type>;

  typename NamedType::Type_ const &value = divisor.value_;
#ifndef NDEBUG
  if (numeric_operations::IsEqual(value, Type{}))
      Throw(DomainError("Division by zero (" + std::to_string(value) + ") encountered."),
            "splinelib::sources::utilities::operator/::NamedType");
#endif
  return NamedType{dividend / value};
}

template<typename Name, typename Type>
constexpr NamedType<Name, Type> operator/(NamedType<Name, Type> const &dividend, Type const &divisor) {
#ifndef NDEBUG
  if (numeric_operations::IsEqual(divisor, Type{}))
      Throw(DomainError("Division by zero (" + std::to_string(divisor) + ") encountered."),
            "splinelib::sources::utilities::operator/::NamedType");
#endif
  return NamedType<Name, Type>{dividend.value_ / divisor};
}

template<typename Name, typename Type>
constexpr bool IsEqual(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs, Type const &tolerance) {
#ifndef NDEBUG
  try {
    numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const &exception) {
    Throw(exception, "splinelib::sources::utilities::IsEqual::NamedType");
  }
#endif
  return numeric_operations::IsEqual(lhs.value_, rhs.value_, tolerance);
}

template<typename Name, typename Type>
constexpr bool IsLessOrEqual(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs, Type const &tolerance)
{
#ifndef NDEBUG
  try {
    numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const &exception) {
    Throw(exception, "splinelib::sources::utilities::IsLessOrEqual::NamedType");
  }
#endif
  return IsEqual(lhs, rhs, tolerance) ? true : numeric_operations::IsLess(lhs.value_, rhs.value_, tolerance);
}

template<typename Name, typename Type>
constexpr bool IsGreaterOrEqual(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs,
                                Type const &tolerance) {
#ifndef NDEBUG
  try {
    numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const &exception) {
    Throw(exception, "splinelib::sources::utilities::IsGreaterOrEqual::NamedType");
  }
#endif
  return !IsLess(lhs, rhs, tolerance);
}

template<typename Name, typename Type>
constexpr bool IsLess(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs, Type const &tolerance) {
#ifndef NDEBUG
  try {
    numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const &exception) {
    Throw(exception, "splinelib::sources::utilities::IsLess::NamedType");
  }
#endif
  return IsGreater(rhs, lhs, tolerance);
}

template<typename Name, typename Type>
constexpr bool IsGreater(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs,
                                Type const &tolerance) {
#ifndef NDEBUG
  try {
    numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const &exception) {
    Throw(exception, "splinelib::sources::utilities::IsGreater::NamedType");
  }
#endif
  return !IsLessOrEqual(lhs, rhs, tolerance);
}

template<typename Name, typename Type>
constexpr bool operator==(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs) {
  return IsEqual(lhs, rhs);
}

template<typename Name, typename Type>
constexpr bool operator!=(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs) {
  return (!IsEqual(lhs, rhs));
}

template<typename Name, typename Type>
constexpr bool operator<=(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs) {
  return IsLessOrEqual(lhs, rhs);
}

template<typename Name, typename Type>
constexpr bool operator>=(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs) {
  return IsGreaterOrEqual(lhs, rhs);
}

template<typename Name, typename Type>
constexpr bool operator<(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs) {
  return IsLess(lhs, rhs);
}

template<typename Name, typename Type>
constexpr bool operator>(NamedType<Name, Type> const &lhs, NamedType<Name, Type> const &rhs) {
  return IsGreater(lhs, rhs);
}

template<typename Name, typename Type>
constexpr Type const & NamedType<Name, Type>::Get() const {
  return value_;
}

template<typename Name, typename Type>
constexpr Type & NamedType<Name, Type>::Get() {
  return value_;
}

template<typename Name, typename Type>
void NamedType<Name, Type>::ForEach(int const &first, int const &size,
                                    std::function<void(NamedType<Name, int> const &)> const &function) {
  static_assert(std::is_same_v<Type_, int>, "NamedType must be of integer type!");

  for (NamedType named_integer{first}; named_integer < NamedType{size}; ++named_integer) function(named_integer);
}

#ifndef NDEBUG
template<typename Name, typename Type>
void NamedType<Name, Type>::ThrowIfNamedIntegerIsOutOfBounds(NamedType<Name, int> const &named_integer,
                                                             int const &maximum_value) {
  static_assert(std::is_same_v<Type_, int>, "NamedType must be of integer type!");

  using std::to_string;

  int const &integer = named_integer.Get();
  if (integer > maximum_value)
      throw OutOfRange("The integer " + to_string(integer) + " with name " + typeid(NamedType).name() + " must not be "
                       "greater than " + to_string(maximum_value) + ".");
}
#endif
