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

#ifndef UTILS_COMPUTATIONAL_DIFFERENTIATION_ALGO_DIFF_TYPE_HPP
#define UTILS_COMPUTATIONAL_DIFFERENTIATION_ALGO_DIFF_TYPE_HPP

// c++ library
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

namespace bezman::utils::computational_differentiation {

/*!
 * @class AlgoDiffType
 *
 * This class implements a simple data type for automatic differentiation. The
 * basic mathematical operations have been overloaded using appropriate rules of
 * differentiation in order to determine the value of a computation as well as
 * the derivative with respect to the input at the same time.
 *
 * This implementation is based on and largely inspired by @danielwolff1's
 * implementation in campiga, modified to be used at runtime. This slows down
 * the code, but provides more flexibility.
 *
 * Supported operations:
 *  - Addition
 *  - Subtraction
 *  - Multiplication
 *  - Division
 *  - log, log10 and exp
 *  - sqrt, power
 *  - abs
 *  - Boolean operations
 *  - sin, cos, tan, asin, acos, atan
 *
 * @tparam Scalar      The intrinsic data type for storing the value and the
 *                     derivative components
 * @tparam numberOfDerivatives
 */
template <typename Scalar, std::size_t numberOfDerivatives = 0>
class AlgoDiffType {
 public:
  /// Alias for the intrinsic scalar data type
  using Scalar_ = Scalar;

  /// Alias for the indexing
  using IndexingType_ = std::size_t;

  /// Alias for the data type of this class
  using ADT_ = AlgoDiffType<Scalar, numberOfDerivatives>;

  /// Alias for the data type storing the gradient
  using DerivType_ =
      std::conditional_t<numberOfDerivatives == 0, std::vector<Scalar_>,
                         std::array<Scalar_, numberOfDerivatives>>;

 private:
  /*!
   * Constructor which initializes all member variables directly.
   * This is only meant to be used by this class internally.
   *
   * @param v    The value
   * @param d    The gradient
   */
  AlgoDiffType(const Scalar_ &v, const DerivType_ &d) : v_(v), d_(d) {}

  /*!
   * Stores the result of the computation carried out with the
   * i-th component of vector \f$ x \f$
   */
  Scalar v_;

  /*!
   * Stores the gradient \f$ \nabla_{x_i} f(x) \f$ of the computation carried
   * out with respect to component \f$ x_i \f$
   */
  DerivType_ d_;

 public:
  /// Default Copy Constructor
  constexpr AlgoDiffType(const ADT_ &rhs) = default;

  /// Default Move Constructor
  constexpr AlgoDiffType(ADT_ &&rhs) = default;

  /// Empty constructor which initializes all member variables to zero
  AlgoDiffType() : v_{}, d_{} {}

  /// Scalar Constructor without Derivative
  template <std::size_t dummy = numberOfDerivatives,
            typename std::enable_if_t<dummy == 0> * = nullptr>
  AlgoDiffType(const Scalar_ &value, const IndexingType_ &n_derivatives)
      : v_{value},
        d_(n_derivatives, Scalar_{}){};  // d_{size} initializes to Scalar{}

  /// Scalar Constructor without Derivative
  template <std::size_t dummy = numberOfDerivatives,
            typename std::enable_if_t<dummy != 0> * = nullptr>
  AlgoDiffType(const Scalar_ &value)
      : v_{value}, d_{} {};  // d_{size} initializes to Scalar{}

  /// Define a destructor
  ~AlgoDiffType() = default;

  /*!
   * Constructor for user initialization of an instance. It initializes an
   * instance with the value \f$ x_i \f$ for which the computation is supposed
   * to be carried out and sets its derivative to the canonical basis vector \f$
   * e_i \f$, which corresponds to the component with respect to which the
   * derivative is supposed to be computed.
   */
  template <std::size_t dummy = numberOfDerivatives,
            typename std::enable_if_t<dummy == 0> * = nullptr>
  AlgoDiffType(const Scalar &value, const IndexingType_ &n_derivatives,
               const IndexingType_ active_component)
      : v_{value}, d_(n_derivatives, Scalar_{}) {
    assert(
        ("Requested derivative out of range.", active_component < d_.size()));
    SetActiveComponent(active_component);
  }

  /*!
   * Constructor for user initialization of an instance. It initializes an
   * instance with the value \f$ x_i \f$ for which the computation is supposed
   * to be carried out and sets its derivative to the canonical basis vector \f$
   * e_i \f$, which corresponds to the component with respect to which the
   * derivative is supposed to be computed.
   */
  template <std::size_t dummy = numberOfDerivatives,
            typename std::enable_if_t<dummy != 0> * = nullptr>
  AlgoDiffType(const Scalar &value, const IndexingType_ active_component)
      : v_{value}, d_{} {
    assert(
        ("Requested derivative out of range.", active_component < d_.size()));
    SetActiveComponent(active_component);
  }

  /*!
   * Overload of assign operator
   *
   * If second derivatives are computed (through nesting AlgoDiffTypes), this
   * overload is necessary to assign ADTs by passing scalars.
   */
  template <typename BaseType>
  AlgoDiffType &operator=(const BaseType &base) {
    *this = AlgoDiffType(base);
    return *this;
  }

  AlgoDiffType &operator=(const AlgoDiffType &t) = default;

  /** @defgroup GetterSetters Getter and Setter Methods
   * Setter methods for values and derivatives
   * @{
   */

  /*!
   * Marks the current variable as the component-th active variable of a vector
   * by first filling the derivative vector with zeros and then setting the
   * provided component to one.
   * @param component    The component in range {0,...,n_derivs-1} which is
   * represented by this variable
   * @note This method is meant for initialization and should not be used within
   *       mathematical computations
   */
  void SetActiveComponent(const size_t component) {
    std::fill(d_.begin(), d_.end(), Scalar_{});
    d_[component] = 1.0;
  }

  /// Getter Function
  constexpr Scalar GetValue() const { return v_; }

  /// Getter derivative
  constexpr DerivType_ GetDerivatives() const { return d_; }

  /*!
   * Returns the spatial dimensionality of vector \f$ x \f$
   * @return     The value of the template parameter n_derivs
   */
  constexpr unsigned int GetNumberOfDerivatives() const { return d_.size(); }

  /** @} */  // End of Getters and Setters

  /** @defgroup Operators Overloaded Basic Operators
   * All basic operations that change the value or the derivative of the object
   * and that can be overloaded as class members
   * @{
   */

  /// Addition with AlgoDiffType
  constexpr ADT_ operator+(const ADT_ &b) const;

  /// Addition and self assignment with AlgoDiffType
  constexpr ADT_ &operator+=(const ADT_ &b);

  /// Addition with Scalar
  constexpr ADT_ operator+(const Scalar &b) const;

  /// Addition and self assignment with Scalar
  constexpr ADT_ &operator+=(const Scalar &b);

  /// Negate the value of AlgoDiffType
  constexpr ADT_ operator-() const;

  /// Subtraction with AlgoDiffType
  constexpr ADT_ operator-(const ADT_ &b) const;

  /// Subtraction and self assignment with AlgoDiffType
  constexpr ADT_ &operator-=(const ADT_ &b);

  /// Subtraction with Scalar
  constexpr ADT_ operator-(const Scalar &b) const;

  /// Subtraction and self assignment with Scalar
  constexpr ADT_ &operator-=(const Scalar &b);

  /// Multiplication with AlgoDiffType
  constexpr ADT_ operator*(const ADT_ &b) const;

  /// Multiplication and self assignment with AlgoDiffType
  constexpr ADT_ &operator*=(const ADT_ &b);

  /// Multiplication with Scalar
  constexpr ADT_ operator*(const Scalar &b) const;

  /// Multiplication and self assignment with Scalar
  constexpr ADT_ &operator*=(const Scalar &b);

  /// Division by AlgoDiffType
  constexpr ADT_ operator/(const ADT_ &b) const;

  /// Division and self assignment by AlgoDiffType
  constexpr ADT_ &operator/=(const ADT_ &b);

  /// Division by Scalar
  constexpr ADT_ operator/(const Scalar &b) const;

  /// Division and self assignment by Scalar
  constexpr ADT_ &operator/=(const Scalar &b);

  /** @} */  // End of Basic operations

  /** @defgroup BoolOperators Boolean Operators
   * All boolean operations that change that can be overloaded as class members
   * @{
   */

  /// Greater operator
  constexpr bool operator>(const Scalar &b) const { return v_ > b; };

  /// Greater operator
  constexpr bool operator>(const ADT_ &b) const { return v_ > b.v_; };

  /// Greater equal operator
  constexpr bool operator>=(const Scalar &b) const { return v_ >= b; };

  /// Greater equal operator
  constexpr bool operator>=(const ADT_ &b) const { return v_ >= b.v_; };

  /// Smaller operator
  constexpr bool operator<(const Scalar &b) const { return v_ < b; };

  /// Smaller operator, delegate to operator>=(const ADT_ &adt)
  constexpr bool operator<(const ADT_ &b) const { return v_ < b; };

  /// Smaller equal operator
  constexpr bool operator<=(const Scalar &b) const { return v_ <= b; };

  /// Smaller equal operator
  constexpr bool operator<=(const ADT_ &b) const { return v_ <= b; };

  /// Equal operator
  constexpr bool operator==(const Scalar &b) const { return v_ == b; };

  /// Equal operator for AlgoDiffType (assuming only value considered)
  constexpr bool operator==(const ADT_ &b) const { return v_ == b.v_; };

  /// Inequality operator
  constexpr bool operator!=(const Scalar &b) const { return v_ != b; };

  /// Inequality operator for AlgoDiffType (assuming only value considered)
  constexpr bool operator!=(const ADT_ &b) const { return v_ != b; };

  /** @} */  // End of Boolean Operators

  /** @defgroup friendInjections Friend Injections
   * All operations were the order prohibits the use of class member functions
   * @{
   */

  /// Simple Default Output stream overload
  template <typename ScalarF, std::size_t numberOfDerivativesF>
  friend std::ostream &operator<<(
      std::ostream &os, const AlgoDiffType<ScalarF, numberOfDerivativesF> &a);

  /// Addition
  template <typename ScalarF, std::size_t numberOfDerivativesF>
  friend constexpr AlgoDiffType<ScalarF, numberOfDerivativesF> operator+(
      const ScalarF &a, const AlgoDiffType<ScalarF, numberOfDerivativesF> &b);

  /// Subtraction
  template <typename ScalarF, std::size_t numberOfDerivativesF>
  friend constexpr AlgoDiffType<ScalarF, numberOfDerivativesF> operator-(
      const ScalarF &a, const AlgoDiffType<ScalarF, numberOfDerivativesF> &b);

  /// Multiplication
  template <typename ScalarF, std::size_t numberOfDerivativesF>
  friend constexpr AlgoDiffType<ScalarF, numberOfDerivativesF> operator*(
      const ScalarF &a, const AlgoDiffType<ScalarF, numberOfDerivativesF> &b);

  /// Division
  template <typename ScalarF, std::size_t numberOfDerivativesF>
  friend constexpr AlgoDiffType<ScalarF, numberOfDerivativesF> operator/(
      const ScalarF &a, const AlgoDiffType<ScalarF, numberOfDerivativesF> &b);

  /// Natural exponent of AlgoDiffType (e.g. \f$ \exp{x_i} \f$)
  template <typename ScalarF, std::size_t numberOfDerivativesF>
  friend constexpr AlgoDiffType<ScalarF, numberOfDerivativesF> exp(
      const AlgoDiffType<ScalarF, numberOfDerivativesF> &exponent);

  /// Absolute Value
  template <typename ScalarF, std::size_t numberOfDerivativesF>
  friend constexpr AlgoDiffType<ScalarF, numberOfDerivativesF> abs(
      const AlgoDiffType<ScalarF, numberOfDerivativesF> &base);

  /// Power of AlgoDiffType (e.g. \f$ (x_i)^a \f$)
  template <typename ScalarF, std::size_t numberOfDerivativesF>
  friend constexpr AlgoDiffType<ScalarF, numberOfDerivativesF> pow(
      const AlgoDiffType<ScalarF, numberOfDerivativesF> &base,
      const ScalarF &power);

  /// Power of AlgoDiffType (using \f$ (x)^y = exp(ln(x)y) \f$)
  template <typename ScalarF, std::size_t numberOfDerivativesF>
  friend constexpr AlgoDiffType<ScalarF, numberOfDerivativesF> pow(
      const AlgoDiffType<ScalarF, numberOfDerivativesF> &base,
      const AlgoDiffType<ScalarF, numberOfDerivativesF> &power);

  /// Square root of AlgoDiffType (e.g. \f$ \sqrt{x_i} \f$)
  template <typename ScalarF, std::size_t numberOfDerivativesF>
  friend constexpr AlgoDiffType<ScalarF, numberOfDerivativesF> sqrt(
      const AlgoDiffType<ScalarF, numberOfDerivativesF> &radicand);

  /// Natural logarithm of AlgoDiffType (e.g. \f$ \log{x_i} \f$ )
  template <typename ScalarF, std::size_t numberOfDerivativesF>
  friend constexpr AlgoDiffType<ScalarF, numberOfDerivativesF> log(
      const AlgoDiffType<ScalarF, numberOfDerivativesF> &xi);

  /// Logarithm to base 10 of AlgoDiffType (e.g. \f$ \log_{10}{x_i} \f$)
  template <typename ScalarF, std::size_t numberOfDerivativesF>
  friend constexpr AlgoDiffType<ScalarF, numberOfDerivativesF> log10(
      const AlgoDiffType<ScalarF, numberOfDerivativesF> &a);

  /// Cosine function of AlgoDiffType (e.g. \f$ \cos{x_i} \f$)
  template <typename ScalarF, std::size_t numberOfDerivativesF>
  friend constexpr AlgoDiffType<ScalarF, numberOfDerivativesF> cos(
      const AlgoDiffType<ScalarF, numberOfDerivativesF> &a);

  /// Sine function of AlgoDiffType (e.g. \f$ \sin{x_i} \f$)
  template <typename ScalarF, std::size_t numberOfDerivativesF>
  friend constexpr AlgoDiffType<ScalarF, numberOfDerivativesF> sin(
      const AlgoDiffType<ScalarF, numberOfDerivativesF> &a);

  /// Tangent function of AlgoDiffType (e.g. \f$ \tan{x_i} \f$)
  template <typename ScalarF, std::size_t numberOfDerivativesF>
  friend constexpr AlgoDiffType<ScalarF, numberOfDerivativesF> tan(
      const AlgoDiffType<ScalarF, numberOfDerivativesF> &a);

  /// Inverse Cosine function of AlgoDiffType (e.g. \f$ \acos{x_i} \f$)
  template <typename ScalarF, std::size_t numberOfDerivativesF>
  friend constexpr AlgoDiffType<ScalarF, numberOfDerivativesF> acos(
      const AlgoDiffType<ScalarF, numberOfDerivativesF> &a);

  /// Inverse Sine function of AlgoDiffType (e.g. \f$ \asin{x_i} \f$)
  template <typename ScalarF, std::size_t numberOfDerivativesF>
  friend constexpr AlgoDiffType<ScalarF, numberOfDerivativesF> asin(
      const AlgoDiffType<ScalarF, numberOfDerivativesF> &a);

  /// Inverse Tangent function of AlgoDiffType (e.g. \f$ \atan{x_i} \f$)
  template <typename ScalarF, std::size_t numberOfDerivativesF>
  friend constexpr AlgoDiffType<ScalarF, numberOfDerivativesF> atan(
      const AlgoDiffType<ScalarF, numberOfDerivativesF> &a);

  /// Greater operator with a Scalar
  template <typename ScalarF, std::size_t numberOfDerivativesF>
  friend constexpr bool operator>(
      const ScalarF &scalar,
      const AlgoDiffType<ScalarF, numberOfDerivativesF> &adt);

  /// Greater equal operator  with a Scalar
  template <typename ScalarF, std::size_t numberOfDerivativesF>
  friend constexpr bool operator>=(
      const ScalarF &scalar,
      const AlgoDiffType<ScalarF, numberOfDerivativesF> &adt);

  /// Smaller operator with a Scalar
  template <typename ScalarF, std::size_t numberOfDerivativesF>
  friend constexpr bool operator<(
      const ScalarF &scalar,
      const AlgoDiffType<ScalarF, numberOfDerivativesF> &adt);

  /// Smaller equal operator  with a Scalar
  template <typename ScalarF, std::size_t numberOfDerivativesF>
  friend constexpr bool operator<=(
      const ScalarF &scalar,
      const AlgoDiffType<ScalarF, numberOfDerivativesF> &adt);

  /// Equal operator with a Scalar
  template <typename ScalarF, std::size_t numberOfDerivativesF>
  friend constexpr bool operator==(
      const ScalarF &scalar,
      const AlgoDiffType<ScalarF, numberOfDerivativesF> &adt);

  /// Unequal operator  with a Scalar
  template <typename ScalarF, std::size_t numberOfDerivativesF>
  friend constexpr bool operator!=(
      const ScalarF &scalar,
      const AlgoDiffType<ScalarF, numberOfDerivativesF> &adt);

  /** @} */  // End of Friend Injections

};  // end class AlgoDiffType

#include "bezman/src/utils/computational_differentiation/algo_diff_type.inc"

}  // namespace bezman::utils::computational_differentiation

#endif  // UTILS_COMPUTATIONAL_DIFFERENTIATION_ALGO_DIFF_TYPE_HPP
