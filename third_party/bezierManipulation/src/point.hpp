#ifndef SRC_POINT_HPP
#define SRC_POINT_HPP

#include <array>
#include <iomanip>
#include <iostream>

namespace beziermanipulation {

/*
 * Point class
 *
 * Class inherites most of its functionality from std::array, however, scalar
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
template <unsigned int spatial_dimension, typename BaseType = double>
class Point : public std::array<BaseType, spatial_dimension> {
 public:
  /// Provide data to external users
  constexpr static unsigned int kSpatialDimension = spatial_dimension;

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

  /// Substraction using RAII
  constexpr Point operator-(const Point& rhs) const {
    Point<spatial_dimension, BaseType> difference{(*this)};
    difference -= rhs;
    return difference;
  }

  /// Substraction
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

  /// Facilitate User Output
  friend std::ostream& operator<<(std::ostream& os, const Point& p) {
    os << "[" << std::setw(p.output_precision + 2)
       << std::setprecision(p.output_precision) << p[0];
    for (size_t i{1}; i < spatial_dimension; ++i) {
      os << ", " << std::setw(p.output_precision + 2)
         << std::setprecision(p.output_precision) << p[i];
    }
    os << "]";
    return os;
  }
};  // namespace std::array<BaseType,spatial_dimension>
}  // namespace beziermanipulation

#endif  // SRC_POINT_HPP
