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

#ifndef SOURCES_UTILITIES_STD_CONTAINER_OPERATIONS_HPP_
#define SOURCES_UTILITIES_STD_CONTAINER_OPERATIONS_HPP_

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iterator>
#include <memory>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <vector>

#include "Sources/Utilities/error_handling.hpp"
#include "Sources/Utilities/named_type.hpp"
#include "Sources/Utilities/numeric_operations.hpp"

namespace splinelib {

// Class Template Argument Deduction (CTAD) for aliases possible since C++20, but not yet fully supported by compilers.
template<typename Type, size_t size>
using Array = std::array<Type, size>;
template<typename Type>
using SharedPointer = std::shared_ptr<Type>;
template<typename ...Types>
using Tuple = std::tuple<Types...>;
template<typename Type>
using Vector = std::vector<Type>;

}  // namespace splinelib

// STD container operations such as 1.) checking container types at compile time (currently std::array and std::vector),
// 2.) checked and unchecked (i.e., faster) access to containers in debug and release mode, respectively, 3.)
// transforming containers storing NamedTypes, 4.) generalized comparisons of contained/pointed at values, and 5.) basic
// arithmetic operations.
//
// Example:
//   using NamedInts = Vector<NamedInt>;
//   constexpr bool const &kTrue = is_vector<NamedInts>;
//   NamedInts named_ints{NamedInt{1}, NamedInt{2}, NamedInt{3}}, more_named_ints(named_ints);
//   int const &two = GetValue(TransformNamedTypes<Array<int, 3>>(named_ints), Index{1});
//   GetValue(named_ints, Index{3}) = NamedInt{};  // out_of_range (debug mode) or undefined behavior (release mode).
//   AddAndAssignToFirst(more_named_ints, named_ints);  // more_named_ints = (2 * named_ints).
//   DoesContainEqualValues(more_named_ints, named_ints);  // Yields false as containers do not store equal values.
//   double const &five = EuclidianDistance(Vector{3.0, 4.0}, Vector{0.0, 0.0});
namespace splinelib::sources::utilities::std_container_operations {

template<typename Type>
struct IsArrayStruct : std::false_type {};
template<typename Type, size_t size>
struct IsArrayStruct<Array<Type, size>> : std::true_type {};
template<typename Type>
constexpr bool const is_array{IsArrayStruct<Type>::value};
template<typename Type>
struct IsVectorStruct : std::false_type {};
template<typename Type>
struct IsVectorStruct<Vector<Type>> : std::true_type {};
template<typename Type>
constexpr bool const is_vector{IsVectorStruct<Type>::value};

// Member functions front() and back() have undefined behavior for empty containers.
template<typename ContainerType>
constexpr typename ContainerType::value_type const & GetFront(ContainerType const &container);
template<typename ContainerType>
constexpr typename ContainerType::value_type const & GetBack(ContainerType const &container);
template<typename ContainerType, typename Name>
constexpr typename ContainerType::value_type const & GetValue(ContainerType const &container,
                                                              NamedType<Name, int> const &index);

template<typename ContainerTypeTo, typename ContainerTypeFrom>
constexpr ContainerTypeTo TransformNamedTypes(ContainerTypeFrom const &from);

// Comparison based on given tolerance.
template<typename ContainerType>
constexpr bool DoesContainEqualValues(ContainerType const &lhs, ContainerType const &rhs,
    typename ContainerType::value_type::Type_ const &tolerance =
        numeric_operations::GetEpsilon<typename ContainerType::value_type::Type_>());
template<typename ContainerType>
constexpr bool DoesContainPointersToEqualValues(ContainerType const &lhs, ContainerType const &rhs,
    typename ContainerType::value_type::element_type::Type_ const &tolerance =
    numeric_operations::GetEpsilon<typename ContainerType::value_type::element_type::Type_>());

template<typename ContainerType>
constexpr ContainerType & AddAndAssignToFirst(ContainerType &lhs, ContainerType const &rhs);
template<typename ContainerType, typename ...ContainerTypes>
constexpr ContainerType & AddAndAssignToFirst(ContainerType &lhs, ContainerType const &rhs,
                                              ContainerTypes const &...further_rhs);
template<typename ContainerType, typename ...ContainerTypes>
constexpr ContainerType Add(ContainerType const &lhs, ContainerTypes const &...rhs);
template<typename ContainerType>
constexpr ContainerType & SubtractAndAssignToFirst(ContainerType &lhs, ContainerType const &rhs);
template<typename ContainerType, typename ...ContainerTypes>
constexpr ContainerType & SubtractAndAssignToFirst(ContainerType &lhs, ContainerType const &rhs,
                                                   ContainerTypes const &...further_rhs);
template<typename ContainerType, typename ...ContainerTypes>
constexpr ContainerType Subtract(ContainerType const &lhs, ContainerTypes const &...rhs);
template<typename ContainerType>
constexpr ContainerType Multiply(ContainerType const &container,
                                 typename ContainerType::value_type::Type_ const &factor);
template<typename ContainerType>
ContainerType & DivideAndAssignToFirst(ContainerType &container,
    typename ContainerType::value_type::Type_ const &divisor,
    typename ContainerType::value_type::Type_ const &tolerance =
        numeric_operations::GetEpsilon<typename ContainerType::value_type::Type_>());
template<typename ContainerType>
constexpr ContainerType Divide(ContainerType const &container, typename ContainerType::value_type::Type_ const &divisor,
                               typename ContainerType::value_type::Type_ const &tolerance =
                                   numeric_operations::GetEpsilon<typename ContainerType::value_type::Type_>());

template<typename ContainerType>
constexpr typename ContainerType::value_type DotProduct(ContainerType const &lhs, ContainerType const &rhs);
template<typename ContainerType>
constexpr typename ContainerType::value_type TwoNorm(ContainerType const &container);
template<typename ContainerType>
constexpr typename ContainerType::value_type EuclidianDistance(ContainerType const &lhs, ContainerType const &rhs);

#ifndef NDEBUG
template<typename ContainerTypeLhs, typename ContainerTypeRhs>
void ThrowIfContainerSizesDiffer(ContainerTypeLhs const &lhs, ContainerTypeRhs const &rhs);
#endif

#include "Sources/Utilities/std_container_operations.inc"

}  // namespace splinelib::sources::utilities::std_container_operations

#endif  // SOURCES_UTILITIES_STD_CONTAINER_OPERATIONS_HPP_
