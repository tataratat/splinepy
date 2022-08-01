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

#ifndef SRC_POINT_HPP
#define SRC_POINT_HPP

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace bezman {

/*
 * Point class
 *
 * Class inherits most of its functionality from std::array, however, scalar
 * multiplication and basic operations were added to facilitate the
 * interpolation of control points in the splines class.
 *
 * Note: Class is minimalistic and extended as required
 *
 * @TODO : Write unit tests for this class
 *
 * @tparam spatial_dimension  Spatial dimension of the point
 * @tparam BaseType           Scalar Type describing the coordinates (to impose
 * CD-Types easier)
 */
template <std::size_t spatial_dimension, typename BaseType = double>
class Point : public std::array<BaseType, spatial_dimension> {
 public:
  /// Provide data to external users
  constexpr static unsigned int kSpatialDimension = spatial_dimension;

  /// Provide Type for external use
  using ScalarType = BaseType;

  /// Output precision
  int output_precision{5};

  /// Use default copy constructor
  constexpr Point(const Point&) = default;

  /// Constructor with zero init from std::array
  constexpr Point() = default;

  /// Variadic aggregate initialization (simplified)
  template <typename... scalars>
  explicit constexpr Point(const scalars&... coords)
      : std::array<BaseType, spatial_dimension>{coords...} {
    static_assert(sizeof...(coords) == spatial_dimension,
                  "You are instantiating a Point object with "
                  "more or less than spatial_dimension components.");
  }

  /// Addition using RAII
  constexpr Point operator+(const Point& rhs) const {
    Point<spatial_dimension, BaseType> sum{(*this)};
    sum += rhs;
    return sum;
  }

  /// Addition
  constexpr Point& operator+=(const Point& rhs) {
    for (unsigned int i = 0; i < spatial_dimension; i++) {
      (*this)[i] = (*this)[i] + rhs[i];
    }
    return (*this);
  }

  /// Subtraction using RAII
  constexpr Point operator-(const Point& rhs) const {
    Point<spatial_dimension, BaseType> difference{(*this)};
    difference -= rhs;
    return difference;
  }

  /// Subtraction
  constexpr Point& operator-=(const Point& rhs) {
    for (unsigned int i = 0; i < spatial_dimension; i++) {
      (*this)[i] = (*this)[i] - rhs[i];
    }
    return (*this);
  }

  /// Inversion
  constexpr Point operator-() const {
    Point inverted_point{(*this)};
    return inverted_point * static_cast<BaseType>(-1);
  }

  /// Multiplication with a scalar
  constexpr Point& operator*=(const BaseType& scale) {
    for (unsigned int i = 0; i < spatial_dimension; ++i) (*this)[i] *= scale;
    return (*this);
  }

  /// Multiplication
  constexpr Point operator*(const BaseType& scale) const {
    Point<spatial_dimension, BaseType> product{(*this)};
    product *= scale;
    return product;
  }

  /// Friend injection for reversed arguments (facilitates readability)
  friend constexpr Point operator*(const double& scale, const Point& point) {
    return point * scale;
  }

  /// Division with a scalar
  constexpr Point& operator/=(const BaseType& scale) {
    const BaseType inv_scale{static_cast<BaseType>(1.) / scale};
    (*this) *= inv_scale;
    return (*this);
  }

  /// Division
  constexpr Point operator/(const BaseType& scale) const {
    Point<spatial_dimension, BaseType> division{(*this)};
    division /= scale;
    return division;
  }

  /// Friend injection for reversed arguments (facilitates readability)
  friend constexpr Point operator/(const double& scale, const Point& point) {
    Point<spatial_dimension, BaseType> division{};
    for (unsigned int i = 0; i < spatial_dimension; ++i)
      division[i] = scale / point[i];
    return division;
  }

  /// Scalar Product with another Point of same size
  template <typename Scalar>
  decltype(Scalar{} * BaseType{}) operator*(
      const Point<spatial_dimension, Scalar>& point) const {
    using ScalarReturn = decltype(Scalar{} * BaseType{});
    ScalarReturn result{};
    for (unsigned int i{}; i < spatial_dimension; i++) {
      result += (*this)[i] * point[i];
    }
    return result;
  }

  /// Calculate norm
  BaseType SquaredEuclidianNorm() const {
    BaseType norm{};
    for (unsigned int i{}; i < spatial_dimension; i++) {
      norm += (*this)[i] * (*this)[i];
    }
    return norm;
  }

  /// Calculate norm
  BaseType EuclidianNorm() const {
    return std::sqrt((*this).SquaredEuclidianNorm());
  }

  /// Facilitate User Output
  friend std::ostream& operator<<(std::ostream& os, const Point& p) {
    os << p.toString();
    return os;
  }

  /// Converts Point into a string object
  std::string toString() const {
    std::ostringstream out{};
    out << "[" << std::setw(output_precision + 2)
        << std::setprecision(output_precision) << (*this)[0];
    for (size_t i{1}; i < spatial_dimension; ++i) {
      out << ", " << std::setw(output_precision + 2)
          << std::setprecision(output_precision) << (*this)[i];
    }
    out << "]";
    return out.str();
  }
};  // namespace std::array<BaseType,spatial_dimension>
}  // namespace bezman

#endif  // SRC_POINT_HPP
