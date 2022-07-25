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

#ifndef SOURCES_UTILITIES_MATH_OPERATIONS_HPP_
#define SOURCES_UTILITIES_MATH_OPERATIONS_HPP_

// Math operations (that are not implemented by the standard library) such as 1.) computing binomial coefficients.
//
// Example:
//   int const &four_choose_2 = ComputeBinomialCoefficient(4, 2);  // The binomial coefficient "4 choose 2" equals 6.
namespace splinelib::sources::utilities::math_operations {

int ComputeBinomialCoefficient(int const &number_of_elements_in_set, int const &number_of_elements_in_subset);

}  // namespace splinelib::sources::utilities::math_operations

#endif  // SOURCES_UTILITIES_MATH_OPERATIONS_HPP_
