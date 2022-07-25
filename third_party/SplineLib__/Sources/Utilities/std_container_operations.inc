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

template<typename ContainerType>
constexpr typename ContainerType::value_type const & GetFront(ContainerType const &container) {
#ifndef NDEBUG
  return container.at(0);
#else
  return container.front();
#endif
}

template<typename ContainerType>
constexpr typename ContainerType::value_type const & GetBack(ContainerType const &container) {
#ifndef NDEBUG
  return container.at(container.size() - 1);
#else
  return container.back();
#endif
}

template<typename ContainerType, typename Name>
constexpr typename ContainerType::value_type const & GetValue(ContainerType const &container,
                                                              NamedType<Name, int> const &index) {
#ifndef NDEBUG
  return container.at(index.Get());
#else
  return container[index.Get()];
#endif
}

template<typename ContainerTypeTo, typename ContainerTypeFrom>
constexpr ContainerTypeTo TransformNamedTypes(ContainerTypeFrom const &from) {
  using TypeFrom = typename ContainerTypeFrom::value_type;
  using TypeTo = typename ContainerTypeTo::value_type;
  using std::transform;

  static_assert(is_named_type<TypeFrom>, "Source standard library container must contain NamedTypes.");

  std::function<TypeTo (TypeFrom const &)> const &transformation =
      [] (TypeFrom const &value) { return TypeTo{static_cast<TypeTo>(value.Get())}; };
  ContainerTypeTo to;
  if constexpr (is_array<ContainerTypeTo>) {
#ifndef NDEBUG
    ThrowIfContainerSizesDiffer(from, to);
#endif
    transform(from.begin(), from.end(), to.begin(), transformation);
  } else if constexpr (is_vector<ContainerTypeTo>) {
    to.reserve(from.size());
    transform(from.begin(), from.end(), std::back_inserter(to), transformation);
  } else {
    static_assert(error_handling::is_false<ContainerTypeTo>, "Only std::array and std::vector are supported.");
  }
  return to;
}

template<typename ContainerType>
constexpr bool DoesContainEqualValues(ContainerType const &lhs, ContainerType const &rhs,
                                      typename ContainerType::value_type::Type_ const &tolerance) {
  using ValueType = typename ContainerType::value_type;
  using std::begin, std::end;

#ifndef NDEBUG
  try {
    numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const &exception) {
    Throw(exception, "splinelib::sources::utilities::std_container_operations::DoesContainEqualValues");
  }
#endif
  return std::equal(begin(lhs), end(lhs), begin(rhs), end(rhs), [&] (ValueType const &lhs_value,
                        ValueType const &rhs_value) { return IsEqual(lhs_value, rhs_value, tolerance); });
}

template<typename ContainerType>
constexpr bool DoesContainPointersToEqualValues(ContainerType const &lhs, ContainerType const &rhs,
    typename ContainerType::value_type::element_type::Type_ const &tolerance) {
  using ValueType = typename ContainerType::value_type;
  using std::begin, std::end;

#ifndef NDEBUG
  try {
    numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const &exception) {
    Throw(exception, "splinelib::sources::utilities::std_container_operations::DoesContainPointersToEqualValues");
  }
#endif
  return std::equal(begin(lhs), end(lhs), begin(rhs), end(rhs),
                    [&] (ValueType const &lhs_pointer, ValueType const &rhs_pointer) { return IsEqual(*lhs_pointer,
                                                                                           *rhs_pointer, tolerance); });
}

template<typename ContainerType>
constexpr ContainerType & AddAndAssignToFirst(ContainerType &lhs, ContainerType const &rhs) {
#ifndef NDEBUG
  if constexpr (is_vector<ContainerType>)
      try {
        ThrowIfContainerSizesDiffer(lhs, rhs);
      } catch (InvalidArgument const &exception) {
        Throw(exception, "splinelib::sources::utilities::std_container_operations::AddAndAssignToFirst");
      }
#endif
  std::transform(lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(), std::plus{});
  return lhs;
}
template<typename ContainerType, typename ...ContainerTypes>
constexpr ContainerType & AddAndAssignToFirst(ContainerType &lhs, ContainerType const &rhs,
                                            ContainerTypes const &...further_rhs) {
  return AddAndAssignToFirst(AddAndAssignToFirst(lhs, rhs), further_rhs...);
}

template<typename ContainerType, typename ...ContainerTypes>
constexpr ContainerType Add(ContainerType const &lhs, ContainerTypes const &...rhs) {
  ContainerType sum(lhs);
#ifndef NDEBUG
  try {
#endif
    AddAndAssignToFirst(sum, rhs...);
#ifndef NDEBUG
  } catch (InvalidArgument const &exception) {
    Throw(exception, "splinelib::sources::utilities::std_container_operations::Add");
  }
#endif
  return sum;
}

template<typename ContainerType>
constexpr ContainerType & SubtractAndAssignToFirst(ContainerType &lhs, ContainerType const &rhs) {
#ifndef NDEBUG
  if constexpr (is_vector<ContainerType>)
      try {
        ThrowIfContainerSizesDiffer(lhs, rhs);
      } catch (InvalidArgument const &exception) {
        Throw(exception, "splinelib::sources::utilities::std_container_operations::SubtractAndAssignToFirst");
      }
#endif
  std::transform(lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(), std::minus{});
  return lhs;
}
template<typename ContainerType, typename ...ContainerTypes>
constexpr ContainerType & SubtractAndAssignToFirst(ContainerType &lhs, ContainerType const &rhs,
                                                   ContainerTypes const &...further_rhs) {
  return SubtractAndAssignToFirst(SubtractAndAssignToFirst(lhs, rhs), further_rhs...);
}

template<typename ContainerType, typename ...ContainerTypes>
constexpr ContainerType Subtract(ContainerType const &lhs, ContainerTypes const &...rhs) {
  ContainerType difference{lhs};
#ifndef NDEBUG
  try {
#endif
    SubtractAndAssignToFirst(difference, rhs...);
#ifndef NDEBUG
  } catch (InvalidArgument const &exception) {
    Throw(exception, "splinelib::sources::utilities::std_container_operations::Subtract");
  }
#endif
  return difference;
}

template<typename ContainerType>
constexpr ContainerType Multiply(ContainerType const &container,
                                 typename ContainerType::value_type::Type_ const &factor) {
  using std::bind, std::multiplies, std::placeholders::_1, std::transform;

  ContainerType multiplied;
  if constexpr (is_array<ContainerType>) {
    transform(container.begin(), container.end(), multiplied.begin(), bind(multiplies{}, _1, factor));
  } else if constexpr (is_vector<ContainerType>) {
    multiplied.reserve(container.size());
    transform(container.begin(), container.end(), std::back_inserter(multiplied), bind(multiplies{}, _1, factor));
  } else {
    static_assert(error_handling::is_false<ContainerType>, "Only std::array and std::vector are supported.");
  }
  return multiplied;
}


template<typename ContainerType>
ContainerType & DivideAndAssignToFirst(ContainerType &container,
                                       typename ContainerType::value_type::Type_ const &divisor,
                                       typename ContainerType::value_type::Type_ const &tolerance) {
#ifndef NDEBUG
  if (numeric_operations::IsEqual(divisor, typename ContainerType::value_type::Type_{}, tolerance))
      Throw(InvalidArgument("Division by zero (" + std::to_string(divisor) + ") encountered."),
            "splinelib::sources::utilities::std_container_operations::DivideAndAssignToFirst");
#endif
  std::transform(container.begin(), container.end(), container.begin(), std::bind(std::divides{}, std::placeholders::_1,
                                                                                  divisor));
  return container;
}

template<typename ContainerType>
constexpr ContainerType Divide(ContainerType const &container, typename ContainerType::value_type::Type_ const &divisor,
                               typename ContainerType::value_type::Type_ const &tolerance) {
  ContainerType divided{container};
#ifndef NDEBUG
  try {
#endif
    DivideAndAssignToFirst(divided, divisor, tolerance);
#ifndef NDEBUG
  } catch (InvalidArgument const &exception) {
    Throw(exception, "splinelib::sources::utilities::std_container_operations::Divide");
  }
#endif
  return divided;
}

template<typename ContainerType>
constexpr typename ContainerType::value_type DotProduct(ContainerType const &lhs, ContainerType const &rhs) {
#ifndef NDEBUG
  if constexpr (is_vector<ContainerType>)
      try {
        ThrowIfContainerSizesDiffer(lhs, rhs);
      } catch (InvalidArgument const &exception) {
        Throw(exception, "splinelib::sources::utilities::std_container_operations::DotProduct");
      }
#endif
  return std::inner_product(lhs.begin(), lhs.end(), rhs.begin(), typename ContainerType::value_type{});
}

template<typename ContainerType>
constexpr typename ContainerType::value_type TwoNorm(ContainerType const &container) {
  return typename ContainerType::value_type{std::sqrt(DotProduct(container, container).Get())};
}

template<typename ContainerType>
constexpr typename ContainerType::value_type EuclidianDistance(ContainerType const &lhs, ContainerType const &rhs) {
#ifndef NDEBUG
  if constexpr (is_vector<ContainerType>)
      try {
        ThrowIfContainerSizesDiffer(lhs, rhs);
      } catch (InvalidArgument const &exception) {
        Throw(exception, "splinelib::sources::utilities::std_container_operations::EuclidianDistance");
      }
#endif
  return TwoNorm(Subtract(lhs, rhs));
}

#ifndef NDEBUG
template<typename ContainerTypeLhs, typename ContainerTypeRhs>
void ThrowIfContainerSizesDiffer(ContainerTypeLhs const &lhs, ContainerTypeRhs const &rhs) {
  using std::to_string;

  int const &size_lhs = lhs.size(), &size_rhs = rhs.size();
  if (size_lhs != size_rhs)
      throw InvalidArgument("Expected equal sizes of the lhs and rhs containers but encountered sizes of " +
                            to_string(size_lhs) + " and " + to_string(size_rhs) + ", respectively.");
}
#endif
