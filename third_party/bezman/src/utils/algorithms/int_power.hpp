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

#ifndef UTILS_ALGORITHMS_INT_POWER_HPP
#define UTILS_ALGORITHMS_INT_POWER_HPP

namespace bezman::utils::algorithms {

/**
 * @brief Recursive function of integer powers
 */
template <typename IntegralType,
          typename E = std::enable_if_t<std::is_integral_v<IntegralType>>>
constexpr IntegralType IntPower(const IntegralType& base,
                                const IntegralType& power) {
  if (power == 0) return 1;
  if (power == 1) return base;

  int tmp = IntPower(base, power / 2);
  if ((power % 2) == 0) {
    return tmp * tmp;
  } else {
    return base * tmp * tmp;
  }
}
}  // namespace bezman::utils::sort

#endif  // UTILS_ALGORITHMS_INT_POWER_HPP
