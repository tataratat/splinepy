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

#ifndef SOURCES_UTILITIES_STRING_OPERATIONS_HPP_
#define SOURCES_UTILITIES_STRING_OPERATIONS_HPP_

#include <algorithm>
// NOLINTNEXTLINE(whitespace/comments)
//#include <format>  // TODO(all): C++20 library providing std::format to replace Write (not yet implemented).
#include <functional>
#include <iterator>
#include <sstream>
#include <string>
#include <type_traits>

#include "Sources/Utilities/error_handling.hpp"
#include "Sources/Utilities/named_type.hpp"
#include "Sources/Utilities/numeric_operations.hpp"
#include "Sources/Utilities/std_container_operations.hpp"

// String operations such as 1.) analyzing whether strings start or end with some string pattern, 2.) manipulating
// strings (appending and trimming), 3.) reading and writing numbers.
//
// Example:
//   StartsWith("hello world", "hello");  // Evaluates to true as "hello world" starts with "hello".
//   String delimited_string{"delimited"};
//   Append(delimited_string, " ", Vector{"by", TrimCharacter(" spaces ", ' ')});
//   StringVector const &splitted_at_spaces = SplitAtDelimiter(delimited_string, ' ');  // {"delimited", "by", "spaces"}
//   Degree const &four = ConvertToNumbers<Degree>("4,3", ',')[0];
//   StringArray const &written = Write<StringArray>(Array{Coordinate{1.8765}}, Precision{4});  // {"1.877"}
namespace splinelib::sources::utilities::string_operations {

template<typename Type>
using Numbers = Vector<Type>;
using String = std::string;
template<int length>
using StringArray = Array<String, length>;
using StringVector = Vector<String>;
using StringVectorConstIterator = StringVector::const_iterator;

bool EndsWith(String const &string, String const &pattern);
bool StartsWith(String const &string, String const &pattern);

template<typename SingleStringOrMultipleStrings>
void Append(String &string, String const &delimiter, SingleStringOrMultipleStrings strings);
String TrimCharacter(String string, char const &character);

StringVector SplitAtDelimiter(String string, char const &delimiter);
template<typename Type>
Type ConvertToNumber(String const &string);
template<typename Type>
Numbers<Type> ConvertToNumbers(String delimited_string, char const &delimiter);
template<typename Type>
String Write(Type const &value,
             Precision const &precision = Precision{numeric_operations::GetPrecision<typename Type::Type_>()});
template<typename ContainerTypeTo, typename ContainerTypeFrom>
ContainerTypeTo Write(ContainerTypeFrom const &from, Precision const &precision =
                          Precision{numeric_operations::GetPrecision<typename ContainerTypeFrom::value_type::Type_>()});

#include "Sources/Utilities/string_operations.inc"

}  // namespace splinelib::sources::utilities::string_operations

namespace splinelib {

using sources::utilities::string_operations::String, sources::utilities::string_operations::StringArray,
      sources::utilities::string_operations::StringVector;

}  // namespace splinelib

#endif  // SOURCES_UTILITIES_STRING_OPERATIONS_HPP_
