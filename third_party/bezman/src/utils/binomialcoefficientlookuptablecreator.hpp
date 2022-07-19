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

#ifndef SRC_UTILS_BINOMIALCOEFFICIENTLOOKUPTABLECREATOR_HPP
#define SRC_UTILS_BINOMIALCOEFFICIENTLOOKUPTABLECREATOR_HPP

#include <array>

namespace bezman::utils {

/*
 * Class to create a look-up table for binomial coefficients
 *
 * @tparam Scalar intrinsic type used to create the lookup table
 * @tparam max_degree maximum depth of the binomial coefficient
 */
template <typename Scalar, unsigned int max_degree>
class BinomialCoefficientLookupTableCreator {
 public:
  static constexpr std::size_t n_values = max_degree * (max_degree + 1) / 2;

  static constexpr std::array<Scalar, n_values> calculate_table() {
    std::array<Scalar, n_values> ret_val{};
    ret_val[0] = static_cast<Scalar>(1);
    for (std::size_t i{0}; i < max_degree; i++) {
      std::size_t offset = i * (i + 1) / 2;
      for (std::size_t j{0}; j <= i; j++) {
        if (i == j || j == 0) {
          ret_val[offset + j] = static_cast<Scalar>(1);
        } else {
          ret_val[offset + j] =
              ret_val[offset - i + j - 1] + ret_val[offset - i + j];
        }
      }
    }
    return ret_val;
  }
};

}

#endif  // SRC_UTILS_BINOMIALCOEFFICIENTLOOKUPTABLECREATOR_HPP
