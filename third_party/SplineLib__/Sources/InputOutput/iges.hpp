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

#ifndef SOURCES_INPUTOUTPUT_IGES_HPP_
#define SOURCES_INPUTOUTPUT_IGES_HPP_

#include "Sources/InputOutput/operations.hpp"
#include "Sources/Utilities/named_type.hpp"
#include "Sources/Utilities/string_operations.hpp"

// Input & output for IGES (cf. IGES 5.3 at, e.g., <http://paulbourke.net/dataformats/iges/IGES.pdf>).
//
// Example:
//   Splines const &splines = Read("input.iges");
namespace splinelib::sources::input_output::iges {

Splines Read(String const &file_name);
void Write(Splines const &splines, String const &file_name, Precision const &precision = kPrecision,
           Tolerance const &tolerance = kEpsilon);

}  // namespace splinelib::sources::input_output::iges

#endif  // SOURCES_INPUTOUTPUT_IGES_HPP_
