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

#ifndef UTILS_TYPE_TRAITS_IS_POINT_HPP
#define UTILS_TYPE_TRAITS_IS_POINT_HPP

#include "bezman/src/point.hpp"

namespace bezman::utils::type_traits
{
  
/// Checker if a template type is an instance of Point
template <typename PointType>
struct isPoint {
  constexpr static bool value = false;
};

/// Checker if a template type is an instance of Point
template <std::size_t spatial_dimension, typename BaseType>
struct isPoint<Point<spatial_dimension, BaseType>> {
  constexpr static bool value = true;
};

/// Alias std naming conform
template <typename PointType>
inline constexpr bool isPoint_v = isPoint<PointType>::value;

} // namespace bezman::utils::type_traits


#endif  // UTILS_TYPE_TRAITS_IS_POINT_HPP
