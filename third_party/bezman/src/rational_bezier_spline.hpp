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

#ifndef SRC_RATIONAL_BEZIER_SPLINE_HPP
#define SRC_RATIONAL_BEZIER_SPLINE_HPP

#include <array>
#include <vector>

#include "bezman/src/bezier_spline.hpp"
#include "bezman/src/point.hpp"

namespace bezman {

// Forward declaration for later use
template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
class RationalBezierSplineGroup;
template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
class BezierSpline;

/*
 * Class describing rational Bezier-Splines
 *
 * @tparam parametric_dimenasion parametric dimension of the spline
 * @tparam PhysicalPointType Type of the control points that are used for
 * interpolations
 * @tparam ScalarType default scalar type used to in the physical domain
 * (e.g. double / AD-Type)
 */
template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType = typename PhysicalPointType::Scalar>
class RationalBezierSpline {
 private:
  // Aliases
  using IndexingType = std::size_t;
  using PointTypePhysical_ = PhysicalPointType;
  using PointTypeParametric_ = Point<parametric_dimension, ScalarType>;
  using PolynomialScalarBezier =
      BezierSpline<parametric_dimension, ScalarType, ScalarType>;
  using PolynomialBezier =
      BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>;

  /*!
   * Numerator Spline, with control points multiplied with their
   * respective weights
   */
  PolynomialBezier weighted_spline_;

  /*!
   * Denominator Spline contains the weight function
   */
  PolynomialScalarBezier weight_function_;

  /**
   * @brief  Wrapper function that checks if the Numerator and Denominator
   * splines are still compatible after operations are performed
   */
  constexpr bool CheckSplineCompatibility() const {
    bool compatible =
        (weighted_spline_.GetDegrees() == weight_function_.GetDegrees());
    return compatible;
  }

  // Friend declarations
  template <std::size_t parent_parametric_dimension,
            typename ParentPhysicalPointType, typename ParentScalarType>
  friend class BezierSpline;
  template <std::size_t parent_parametric_dimension,
            typename ParentPhysicalPointType, typename ParentScalarType>
  friend class RationalBezierSpline;

 public:
  /// Make ScalarType publicly available
  using ScalarType_ = ScalarType;

  /// Make Parametric dimension publicly available
  static constexpr IndexingType kParametricDimensions = parametric_dimension;

  /// Copy constructor
  constexpr RationalBezierSpline(const RationalBezierSpline& bezier_spline) =
      default;

  /// Transform Bezier into rational Bezier
  constexpr RationalBezierSpline(
      const BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
          bezier_spline)
      : weighted_spline_{bezier_spline} {
    // Initialize new weight function spline
    weight_function_ =
        BezierSpline<parametric_dimension, ScalarType, ScalarType>{
            // Degrees can be taken from weighted spline
            weighted_spline_.GetDegrees(),
            // Generate control weights (just 1 for every point)
            std::vector<ScalarType>(
                // Size of control-point vector
                weighted_spline_.NumberOfControlPoints,
                // Value
                static_cast<ScalarType>(1.))};
    assert(CheckSplineCompatibility());
  }

  /// Constructor with two spline - weight and weighted spline
  constexpr RationalBezierSpline(
      const PolynomialBezier& r_weighted_spline_,
      const PolynomialScalarBezier& r_weight_function_)
      : weighted_spline_{r_weighted_spline_},
        weight_function_{r_weight_function_} {
    assert(CheckSplineCompatibility());
  }

  /// Empty constructor
  constexpr RationalBezierSpline() = default;

  /// Constructor with weights
  constexpr RationalBezierSpline(
      const std::array<std::size_t, parametric_dimension> deg)
      : weighted_spline_{deg}, weight_function_{deg} {};

  /// Constructor with control point list
  constexpr RationalBezierSpline(
      const std::array<std::size_t, parametric_dimension> deg,
      const std::vector<PointTypePhysical_> control_point_vector,
      const std::vector<ScalarType> weights)
      : weighted_spline_{deg}, weight_function_{deg, weights} {
    assert(control_point_vector.size() == weights.size());
    assert(weighted_spline.kNumberOfControlPoints == weights.size());
    for (std::size_t i_point{0};
         i_point < weighted_spline_.control_points.size(); i_point++) {
      weighted_spline_.control_points[i_point] =
          control_point_vector[i_point] * weights[i_point];
    }
  }

  /// Move operator
  constexpr RationalBezierSpline& operator=(RationalBezierSpline&& rhs) =
      default;

  /// Move operator
  constexpr RationalBezierSpline& operator=(const RationalBezierSpline& rhs) =
      default;

  /// Getter for Degrees
  constexpr const std::array<std::size_t, parametric_dimension>& GetDegrees()
      const {
    return weighted_spline_.GetDegrees();
  }

  /// Access Control Point Vector directly
  constexpr const std::vector<PhysicalPointType>& GetWeightedControlPoints()
      const {
    return weighted_spline_.control_points;
  }

  /// Access Control Point Vector directly
  constexpr std::vector<PhysicalPointType>& GetWeightedControlPoints() {
    return weighted_spline_.control_points;
  }

  /// Access Weights Vector directly
  constexpr const std::vector<ScalarType>& GetWeights() const {
    return weight_function_.control_points;
  }

  /// Access Weights Vector directly
  constexpr std::vector<ScalarType>& GetWeights() {
    return weight_function_.control_points;
  }

  /// Set Degrees
  constexpr void UpdateDegrees(
      const std::array<std::size_t, parametric_dimension>& new_degrees) {
    weight_function_.UpdateDegrees(new_degrees);
    weighted_spline_.UpdateDegrees(new_degrees);
    return;
  }

  // @NOTE These functions all pass arguments and perform very limited
  // operations

  // Get the weight function spline as const reference
  constexpr const PolynomialScalarBezier& GetWeightFunctionSpline() const {
    return weight_function_;
  }

  /// Retrieve single control point from local indices
  template <typename... T>
  constexpr PointTypePhysical_ ControlPoint(const T... index) const;

  /// Retrieve single control point from local indices (as array)
  constexpr PointTypePhysical_ ControlPoint(
      const std::array<IndexingType, parametric_dimension>& index) const;

  /// Retrieve single weighted eighted control point from local indices
  template <typename... T>
  constexpr const PointTypePhysical_& WeightedControlPoint(
      const T... index) const;

  /// Retrieve single weighted control point from local indices
  template <typename... T>
  constexpr PointTypePhysical_& WeightedControlPoint(const T... index);

  /// Retrieve single weighted control point from local indices (as array)
  constexpr const PointTypePhysical_& WeightedControlPoint(
      const std::array<IndexingType, parametric_dimension>& index) const;

  /// Retrieve single weighted control point from local indices (as array)
  constexpr PointTypePhysical_& WeightedControlPoint(
      const std::array<IndexingType, parametric_dimension>& index);

  /// Retrieve single weight from local indices
  template <typename... T>
  constexpr const ScalarType& Weight(const T... index) const;

  /// Retrieve single weight from local indices
  template <typename... T>
  constexpr ScalarType& Weight(const T... index);

  /// Retrieve single weight from local indices (as array)
  constexpr const ScalarType& Weight(
      const std::array<IndexingType, parametric_dimension>& index) const;

  /// Retrieve single weight from local indices (as array)
  constexpr ScalarType& Weight(
      const std::array<IndexingType, parametric_dimension>& index);

  //------------------- Essential Operations

  /// Evaluate the spline via the deCasteljau algorithm
  constexpr PointTypePhysical_ Evaluate(
      const PointTypeParametric_& par_coords) const;

  /// Evaluate the spline using the de Casteljau algorithm
  template <typename... T>
  constexpr PointTypePhysical_ Evaluate(const T&... par_coords) const {
    return Evaluate(PointTypeParametric_{par_coords...});
  }

  /// Evaluate Basis Functions
  template <typename... T>
  constexpr std::array<std::vector<ScalarType>, parametric_dimension>
  BasisFunctions(const T&... par_coords) const {
    return BasisFunctions(PointTypeParametric_{par_coords...});
  }

  /// Evaluate Basis Functions
  constexpr std::array<std::vector<ScalarType>, parametric_dimension>
  BasisFunctions(const PointTypeParametric_& par_coords) const;

  /// Evaluate Basis Functions
  template <typename... T>
  constexpr std::array<std::vector<ScalarType>, parametric_dimension>
  PolynomialBasisFunctions(const T&... par_coords) const {
    return PolynomialBasisFunctions(PointTypeParametric_{par_coords...});
  }

  /// Evaluate Basis Functions without respecting weights
  constexpr std::array<std::vector<ScalarType>, parametric_dimension>
  PolynomialBasisFunctions(const PointTypeParametric_& par_coords) const;

  /// Evaluate the spline via the explicit precomputation of bernstein values
  constexpr PointTypePhysical_ ForwardEvaluate(
      const PointTypeParametric_& par_coords) const;

  /// Evaluate the spline using explicit precomputation of bernstein values
  template <typename... T>
  constexpr PointTypePhysical_ ForwardEvaluate(const T&... par_coords) const {
    return ForwardEvaluate(PointTypeParametric_{par_coords...});
  }

  //-------------------  Refinement Strategies

  /// Order elevation along a specific parametric dimension
  constexpr RationalBezierSpline& OrderElevateAlongParametricDimension(
      const IndexingType par_dim);

  /*
   * Split the Bezier Spline into two distinct subdivisions
   *
   * Splits the Spline along a specific dimension and returns a group
   * representing the same domain over two splines.
   *
   * @attention Not Implemented
   */
  constexpr RationalBezierSplineGroup<parametric_dimension, PhysicalPointType,
                                      ScalarType>
  SplitAtPosition(const ScalarType& splitting_plane,
                  const IndexingType splitting_dimension = 0) const;

  /**
   * @brief  Derivative along a specific parametric dimension
   *
   * This function is implemented using the Quotient-rule for derivatives
   */
  constexpr RationalBezierSpline DerivativeWRTParametricDimension(
      const IndexingType par_dim) const;

  //-------------------  Operator overloads

  /// Multiplication with scalar
  constexpr RationalBezierSpline& operator*=(const ScalarType& scalar);

  /// Multiplication with scalar RAII
  constexpr RationalBezierSpline operator*(const ScalarType& scalar) const;

  /// Friend injection for reversed order
  friend constexpr RationalBezierSpline operator*(
      const ScalarType& scalar, const RationalBezierSpline& b) {
    return b * scalar;
  }

  /// Addition
  constexpr RationalBezierSpline& operator+=(
      const PhysicalPointType& point_shift);

  /// Addition
  constexpr RationalBezierSpline operator+(
      const PhysicalPointType& point_shift) const;

  /// Friend injection Addition
  friend constexpr RationalBezierSpline operator+(
      const PhysicalPointType& point_shift,
      const RationalBezierSpline& original_spline) {
    return original_spline + point_shift;
  }

  /// Substraction
  constexpr RationalBezierSpline operator-(
      const PhysicalPointType& point_shift) const;

  /// Substraction
  constexpr RationalBezierSpline& operator-=(
      const PhysicalPointType& point_shift);

  /// Inversion
  constexpr RationalBezierSpline operator-() const;

  /// Addition of Two Splines resulting in a new spline that describes the
  /// pointwise addition of the two Beziers
  template <typename PointTypeRHS, typename ScalarRHS>
  constexpr RationalBezierSpline<parametric_dimension,
                                 decltype(PhysicalPointType{} + PointTypeRHS{}),
                                 decltype(ScalarType_{} * ScalarRHS{})>
  operator+(const RationalBezierSpline<parametric_dimension, PointTypeRHS,
                                       ScalarRHS>& rhs) const;

  /// Add two splines of same type
  constexpr RationalBezierSpline& operator+=(const RationalBezierSpline& rhs);

  /// Substraction of Two Splines resulting in a new spline that describes the
  /// pointwise addition of the two Beziers
  template <typename PointTypeRHS, typename ScalarRHS>
  constexpr RationalBezierSpline<parametric_dimension,
                                 decltype(PhysicalPointType{} - PointTypeRHS{}),
                                 decltype(ScalarType_{} * ScalarRHS{})>
  operator-(const RationalBezierSpline<parametric_dimension, PointTypeRHS,
                                       ScalarRHS>& rhs) const;

  /// Subtract two splines of same type
  constexpr RationalBezierSpline& operator-=(const RationalBezierSpline& rhs);

  /// Multiplication of two splines similar to pointwise product
  template <typename PointTypeRHS, typename ScalarRHS>
  constexpr RationalBezierSpline<parametric_dimension,
                                 decltype(PhysicalPointType{} * PointTypeRHS{}),
                                 decltype(ScalarType_{} * ScalarRHS{})>
  operator*(const RationalBezierSpline<parametric_dimension, PointTypeRHS,
                                       ScalarRHS>& rhs) const;

  //------------------- IRIT-related operations

  /// Extract single coordinate spline
  constexpr BezierSpline<parametric_dimension, ScalarType, ScalarType>
  ExtractDimension(const std::size_t& dimension) const;

  /*!
   * Functional Composition between two splines
   *
   * Compose two splines, taking the (*this) spline as the outer funtion and
   * the function argument as the inner function. This works so long as the
   * parametric dimension of the outer function matches the physical dimension
   * of the inner function.
   */
  template <std::size_t parametric_dimension_inner_spline,
            typename PointTypeRHS, typename ScalarRHS>
  constexpr RationalBezierSpline<parametric_dimension_inner_spline,
                                 PhysicalPointType,
                                 decltype(ScalarType_{} * ScalarRHS{})>
  Compose(const BezierSpline<parametric_dimension_inner_spline, PointTypeRHS,
                             ScalarRHS>& inner_function) const {
    return RationalBezierSpline<parametric_dimension_inner_spline,
                                PhysicalPointType,
                                decltype(ScalarType_{} * ScalarRHS{})>{
        // Numerator composition
        weighted_spline_.Compose(inner_function),
        // Denominator composition
        weight_function_.Compose(inner_function)};
  }

  /*!
   * Functional Composition between two rational
   *
   * Compose two splines, taking the (*this) spline as the outer funtion and
   * the function argument as the inner function. This works so long as the
   * parametric dimension of the outer function matches the physical dimension
   * of the inner function.
   *
   * For two rational splines, the weight function cancels out, see
   * BezierSpline::ComposeDenominator
   */
  template <std::size_t parametric_dimension_inner_spline,
            typename PointTypeRHS, typename ScalarRHS>
  constexpr RationalBezierSpline<parametric_dimension_inner_spline,
                                 PhysicalPointType,
                                 decltype(ScalarType_{} * ScalarRHS{})>
  Compose(const RationalBezierSpline<parametric_dimension_inner_spline,
                                     PointTypeRHS, ScalarRHS>& inner_function)
      const {
    return RationalBezierSpline<parametric_dimension_inner_spline,
                                PhysicalPointType,
                                decltype(ScalarType_{} * ScalarRHS{})>{
        // Numerator composition
        weighted_spline_.ComposeNumeratorSpline(inner_function),
        // Denominator composition
        weight_function_.ComposeNumeratorSpline(inner_function)};
  }

  /// Check if can be used for composition
  constexpr bool FitsIntoUnitCube() const {
    return weighted_spline_.FitsIntoUnitCube();
  };

  /// Get maximum restricting corner of spline
  constexpr PhysicalPointType MaximumCorner() const {
    return weighted_spline_.MaximumCorner();
  };

  /// Get minimum restricting corner of spline
  constexpr PhysicalPointType MinimumCorner() const {
    return weighted_spline_.MinimumCorner();
  };
};

#include "bezman/src/rational_bezier_spline.inc"

}  // namespace bezman

#endif  // SRC_RATIONAL_BEZIER_SPLINE_HPP
