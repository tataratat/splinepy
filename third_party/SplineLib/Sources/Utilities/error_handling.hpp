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

#ifndef SOURCES_UTILITIES_ERROR_HANDLING_HPP_
#define SOURCES_UTILITIES_ERROR_HANDLING_HPP_

#include <stdexcept>
#include <string>
#include <type_traits>

// Error handling such as 1.) tracing exceptions when compiled in debug mode and 2.) providing a dependent false for
// static_asserts.
//
// Example:
//   try {
//     FunctionThrowingRuntimeError();  // Throwing function call within function named CallerOfThrowingFunction.
//   } catch (RuntimeError const &exception) {
//     Throw(exception, "splinelib::sources::utilities::CallerOfThrowingFunction");  // Propagating the exception.
//   }
namespace splinelib::sources::utilities::error_handling {

template<typename Type>
struct DependentFalse : std::false_type {};
template<typename Type>
constexpr bool const is_false{DependentFalse<Type>::value};

#ifndef NDEBUG
using DomainError = std::domain_error;
using InvalidArgument = std::invalid_argument;
using Message = std::string;
using OutOfRange = std::out_of_range;
using RuntimeError = std::runtime_error;

template<typename Exception>
void Throw(Exception const &exception, Message const &name, int const &dimension = -1);

#include "Sources/Utilities/error_handling.inc"
#endif

}  // namespace splinelib::sources::utilities::error_handling

#ifndef NDEBUG
namespace splinelib {

using sources::utilities::error_handling::DomainError, sources::utilities::error_handling::InvalidArgument,
      sources::utilities::error_handling::Message, sources::utilities::error_handling::OutOfRange,
      sources::utilities::error_handling::RuntimeError;
using sources::utilities::error_handling::Throw;

}  // namespace splinelib
#endif

#endif  // SOURCES_UTILITIES_ERROR_HANDLING_HPP_
