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

#ifndef SOURCES_INPUTOUTPUT_VTK_HPP_
#define SOURCES_INPUTOUTPUT_VTK_HPP_

#include "Sources/InputOutput/operations.hpp"
#include "Sources/Utilities/named_type.hpp"
#include "Sources/Utilities/std_container_operations.hpp"
#include "Sources/Utilities/string_operations.hpp"

// Output for VTK (cf. legacy format at <https://www.kitware.com/products/books/VTKUsersGuide.pdf>).
//
// Example:
//   Write({curve}, "out.vtk", {{NumbersOfParametricCoordinates::value_type::value_type{10}}});  // 9 elements.
namespace splinelib::sources::input_output::vtk {

using NumbersOfParametricCoordinates = Vector<Vector<Length>>;

void Sample(Splines const &splines, String const &file_name,
            NumbersOfParametricCoordinates const &numbers_of_parametric_coordinates,
            Tolerance const &tolerance = kEpsilon, Precision const &precision = kPrecision);

}  // namespace splinelib::sources::input_output::vtk

#endif  // SOURCES_INPUTOUTPUT_VTK_HPP_
