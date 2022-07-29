/*
MIT License

Copyright (c) 2022 zwar@ilsb.tuwien.ac.at

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef SRC_UTILS_FASTBINOMIALCOEFFICIENT_HPP
#define SRC_UTILS_FASTBINOMIALCOEFFICIENT_HPP

#include <array>

#include "bezman/src/utils/binomialcoefficientlookuptablecreator.hpp"

#ifndef MAX_BINOMIAL_DEGREE
#define MAX_BINOMIAL_DEGREE 30u
#endif

namespace bezman::utils {

/*
 * Lookup Class with static binomial coeffient
 *
 * Create a Lookup Table at compile time and provides fast and efficient look up
 * of binomial coefficient. Maximum degree directly correlates with the maximum
 * degree of the bezier spline to be considered.
 */
class FastBinomialCoefficient {
 private:
  /// Sclar Default Type
  using ScalarType = std::size_t;

  /// Size of Lookup table if to large, needs to be removed from stack -> heap
  constexpr static auto n_values =
      BinomialCoefficientLookupTableCreator<ScalarType,
                                            MAX_BINOMIAL_DEGREE>::n_values;

  /// Lookup Table from generator
  constexpr static std::array<ScalarType, n_values> look_up =
      BinomialCoefficientLookupTableCreator<
          ScalarType, MAX_BINOMIAL_DEGREE>::calculate_table();

 public:
  /// Binomial coefficient (careful with asserts in -Ofast mode)
  constexpr static ScalarType choose(const unsigned int n,
                                     const unsigned int i) {
    assert(n < MAX_BINOMIAL_DEGREE);
    return look_up[n * (n + 1) / 2 + i];
  }
};
}  // namespace bezman::utils

#endif  // SRC_UTILS_FASTBINOMIALCOEFFICIENT_HPP
