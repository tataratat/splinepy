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

#include "Sources/Utilities/math_operations.hpp"

//#include <cmath> // NOLINT(whitespace/comments)

#include "Sources/Utilities/error_handling.hpp"

namespace splinelib::sources::utilities::math_operations {

int ComputeBinomialCoefficient(int const &number_of_elements_in_set, int const &number_of_elements_in_subset) {
#ifndef NDEBUG
  using std::to_string;

  Message const kName{"splinelib::sources::utilities::math_operations::ComputeBinomialCoefficient"};

  if (number_of_elements_in_set < 0) {
    Throw(InvalidArgument("The number of elements in the set (" + to_string(number_of_elements_in_set) + ") must not "
                          "be negative."), kName);
  } else if (number_of_elements_in_subset < 0) {
    Throw(InvalidArgument("The number of elements in the subset (" + to_string(number_of_elements_in_subset) + ") must "
                          "not be negative."), kName);
  } else if (number_of_elements_in_subset > number_of_elements_in_set) {
    Throw(InvalidArgument("The number of elements in the subset (" + to_string(number_of_elements_in_subset) + ") must "
                          "not be greater than the number of elements in the set (" +
                          to_string(number_of_elements_in_set) + ")."), kName);
  }
#endif
//  NOLINTNEXTLINE(whitespace/todo)
//  TODO(all): C++17 library cmath should provide std::beta (not yet implemented).
//  int const &number_of_elements_in_set_plus_1 = (number_of_elements_in_set + 1);
//  return (1 / (number_of_elements_in_set_plus_1 *
//      std::beta(number_of_elements_in_set_plus_1 - number_of_elements_in_subset, number_of_elements_in_subset + 1)));
  if ((number_of_elements_in_subset == number_of_elements_in_set) || (number_of_elements_in_subset == 0)) {
    return 1;
  } else {
    int const &number_of_elements_in_set_minus_1 = (number_of_elements_in_set - 1);
    return (ComputeBinomialCoefficient(number_of_elements_in_set_minus_1, number_of_elements_in_subset) +
            ComputeBinomialCoefficient(number_of_elements_in_set_minus_1, number_of_elements_in_subset - 1));
  }
}

}  // namespace splinelib::sources::utilities::math_operations
