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

#ifndef SOURCES_UTILITIES_SYSTEM_OPERATIONS_HPP_
#define SOURCES_UTILITIES_SYSTEM_OPERATIONS_HPP_

#include <ctime>
#include <fstream>
#include <ios>
#include <string>

#include "Sources/Utilities/error_handling.hpp"

// System operations such as 1.) getting the local time and 2.) opening files.
//
// Example:
//   LocalTime const &local_time = GetLocalTime();
//   OutputStream output_stream{Open<OutputStream, kModeOut>("file.out")};
namespace splinelib::sources::utilities::system_operations {

using File = std::string;
using InputStream = std::ifstream;
using LocalTime = std::tm;
using Mode = std::ios_base::openmode;
using OutputStream = std::ofstream;

constexpr sources::utilities::system_operations::Mode const kModeAppend{std::ios_base::app}, kModeIn{std::ios_base::in},
                                                            kModeOut{std::ios_base::out};

LocalTime GetLocalTime();

template<typename FileStream, Mode mode>
FileStream Open(File const &file);

#include "Sources/Utilities/system_operations.inc"

}  // namespace splinelib::sources::utilities::system_operations

#endif  // SOURCES_UTILITIES_SYSTEM_OPERATIONS_HPP_
