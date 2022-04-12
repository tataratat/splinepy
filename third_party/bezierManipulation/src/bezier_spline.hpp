#ifndef SRC_BEZIER_SPLINE_HPP
#define SRC_BEZIER_SPLINE_HPP

#include <algorithm>
#include <array>
#include <cassert>
#include <numeric>
#include <vector>

#include "bezierManipulation/src/point.hpp"
#include "bezierManipulation/src/utils/fastbinomialcoefficient.hpp"

namespace beziermanipulation {

// Forward declaration for later use
template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
class BezierSplineGroup;

/*
 * Class describing BezierSplines
 *
 * Core of the library
 *
 * @tparam parametric_dimenasion parametric dimension of the spline
 * @tparam PhysicalPointType Type of the control points that are used for
 * interpolations
 * @tparam ScalarType default scalar type used to in the physical domain
 * (e.g. double / AD-Type)
 */
template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType = typename PhysicalPointType::Scalar>
class BezierSpline {
 private:
  using IndexingType = std::size_t;
  using PointTypePhysical_ = PhysicalPointType;
  using PointTypeParametric_ = Point<parametric_dimension, ScalarType>;

  constexpr void UpdateIndexOffsets_();

  template <typename... Indices>
  constexpr PhysicalPointType AddUpContributionsToControlPointVector_(
      PhysicalPointType& evaluation_point,
      const std::array<std::array<ScalarType, MAX_BINOMIAL_DEGREE + 1>,
                       parametric_dimension>& factors,
      const ScalarType& factor_product, const Indices&... indices) const;

  /// Polynomial degrees
  std::array<IndexingType, parametric_dimension> degrees{};

 public:
  using ScalarType_ = ScalarType;

  /// Offsets in Row based control point storage
  std::array<IndexingType, parametric_dimension> index_offsets{};
  /// Number of control points
  IndexingType NumberOfControlPoints{};
  /// List of all control points in "Row-based" order
  std::vector<PointTypePhysical_> control_points{};

  /// Make Parametric dimension publicly available
  static constexpr IndexingType kParametricDimensions = parametric_dimension;

  /*
   * Retrieve individual indices
   *
   * @param local_indices single value index as control point index
   */
  constexpr std::array<IndexingType, parametric_dimension> LocalToGlobalIndex(
      const IndexingType& local_index) const;

  /*
   * Calculate the new control point that result from the multiplication
   * between two bezier splines
   * Note : Using the notation in Gershons diss (eq. 2.13)
   */
  template <typename PhysicalPointLHS, typename ScalarLHS,
            typename PhysicalPointRHS, typename ScalarRHS, typename... T>
  constexpr void CombineControlPointsForProduct_(
      const BezierSpline<parametric_dimension, PhysicalPointLHS, ScalarLHS>&
          P_spline,
      const BezierSpline<parametric_dimension, PhysicalPointRHS, ScalarRHS>&
          Q_spline,
      const std::array<IndexingType, parametric_dimension>& ctpsIndex,
      const ScalarType factor, const T&... indices);

  /// Copy constructor
  constexpr BezierSpline(const BezierSpline& bezier_spline) = default;

  /// Empty constructor
  constexpr BezierSpline() = default;

  /// Empty constructor with degrees
  constexpr BezierSpline(
      const std::array<std::size_t, parametric_dimension> deg)
      : degrees{deg} {
    NumberOfControlPoints = 1u;
    for (unsigned int i{}; i < parametric_dimension; i++)
      NumberOfControlPoints *= degrees[i] + 1;
    control_points.resize(NumberOfControlPoints);
    UpdateIndexOffsets_();
  };

  /// Constructor with control point list
  constexpr BezierSpline(
      const std::array<std::size_t, parametric_dimension> deg,
      const std::vector<PointTypePhysical_> points)
      : degrees{deg}, control_points{points} {
    NumberOfControlPoints = 1u;
    for (unsigned int i{}; i < parametric_dimension; i++)
      NumberOfControlPoints *= degrees[i] + 1;

    UpdateIndexOffsets_();
    assert(NumberOfControlPoints == points.size());
  };

  /// Move operator
  constexpr BezierSpline& operator=(BezierSpline&& rhs) = default;

  /// Move operator
  constexpr BezierSpline& operator=(const BezierSpline& rhs) = default;

  /// Getter for Degrees
  constexpr const std::array<std::size_t, parametric_dimension>& GetDegrees()
      const {
    return degrees;
  }

  /// Getter for Degrees
  constexpr std::array<std::size_t, parametric_dimension> GetDegrees() {
    return degrees;
  }

  /// Set Degrees
  constexpr void UpdateDegrees(const std::array<std::size_t, parametric_dimension>& new_degrees){
    degrees = new_degrees;
    for (unsigned int i{}; i < parametric_dimension; i++)
      NumberOfControlPoints *= degrees[i] + 1;
    control_points.resize(NumberOfControlPoints);
    UpdateIndexOffsets_();
  }

  /// Retrieve single control point from local indices
  template <typename... T>
  constexpr const PointTypePhysical_& ControlPoint(const T... index) const;

  /// Retrieve single control point from local indices
  template <typename... T>
  constexpr PointTypePhysical_& ControlPoint(const T... index);

  /// Retrieve single control point from local indices (as array)
  constexpr const PointTypePhysical_& ControlPoint(
      const std::array<IndexingType, parametric_dimension>& index) const;

  /// Retrieve single control point from local indices (as array)
  constexpr PointTypePhysical_& ControlPoint(
      const std::array<IndexingType, parametric_dimension>& index);

  /// Order elevation along a specific parametric dimension
  constexpr BezierSpline& OrderElevateAlongParametricDimension(
      const IndexingType par_dim);

  /// Derivative along a specific parametric dimension
  constexpr BezierSpline DerivativeWRTParametricDimension(
      const IndexingType par_dim) const;

  /// Evaluate the spline using the de Casteljau algorithm
  template <typename... T>
  constexpr PointTypePhysical_ Evaluate(const T&... par_coords) const {
    return Evaluate(PointTypeParametric_{par_coords...});
  }

  /// Evaluate the spline via the deCasteljau algorithm
  constexpr PointTypePhysical_ Evaluate(
      const PointTypeParametric_& par_coords) const;

  /// Evaluate the spline via the explicit precomputation of bernstein values
  constexpr PointTypePhysical_ ForwardEvaluate(
      const PointTypeParametric_& par_coords) const;

  /// Evaluate the spline using explicit precomputation of bernstein values
  template <typename... T>
  constexpr PointTypePhysical_ ForwardEvaluate(const T&... par_coords) const {
    return ForwardEvaluate(PointTypeParametric_{par_coords...});
  }

  /// Addition of Two Splines resulting in a new spline that describes the
  /// pointwise addition of the two Beziers
  template <typename PointTypeRHS, typename ScalarRHS>
  constexpr BezierSpline<parametric_dimension,
                         decltype(PhysicalPointType{} + PointTypeRHS{}),
                         decltype(ScalarType_{} + ScalarRHS{})>
  operator+(const BezierSpline<parametric_dimension, PointTypeRHS, ScalarRHS>&
                rhs) const;

  /// Add two splines of same type
  constexpr BezierSpline& operator+=(BezierSpline rhs);

  /// Check if two splines are equivalent
  constexpr bool operator==(const BezierSpline& rhs) const;

  /// Multiplication of two splines similar to pointwise product
  template <typename PointTypeRHS, typename ScalarRHS>
  constexpr BezierSpline<parametric_dimension,
                         decltype(PhysicalPointType{} * PointTypeRHS{}),
                         decltype(ScalarType_{} * ScalarRHS{})>
  operator*(const BezierSpline<parametric_dimension, PointTypeRHS, ScalarRHS>&
                rhs) const;

  /// Extract single coordinate spline
  constexpr BezierSpline<parametric_dimension, ScalarType, ScalarType>
  ExtractDimension(unsigned int dimension) const;

  constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>
  RaisePower(const unsigned int power) const;

  /// Multiplication with scalar
  constexpr BezierSpline& operator*=(const ScalarType& scalar);

  /// Multiplication with scalar RAII
  constexpr BezierSpline operator*(const ScalarType& scalar) const;

  /// Friend injection for reversed order
  friend constexpr BezierSpline operator*(const ScalarType& scalar,
                                          const BezierSpline& b) {
    return b * scalar;
  }

  /// Check if can be used for composition
  constexpr bool FitsIntoUnitCube() const;

  /// Addition
  constexpr BezierSpline& operator+=(const PhysicalPointType& point_shift);

  /// Addition
  constexpr BezierSpline operator+(const PhysicalPointType& point_shift) const;

  /// Friend injection Addition
  friend constexpr BezierSpline operator+(const PhysicalPointType& point_shift,
                                          const BezierSpline& original_spline) {
    return original_spline + point_shift;
  }

  /// Substraction
  constexpr BezierSpline operator-(const PhysicalPointType& point_shift) const;

  /// Substraction
  constexpr BezierSpline& operator-=(const PhysicalPointType& point_shift);

  /// Inversion
  constexpr BezierSpline operator-() const;

  /// Get maximum restricting corner of spline
  constexpr PhysicalPointType MaximumCorner() const;

  /// Get minimum restricting corner of spline
  constexpr PhysicalPointType MinimumCorner() const;

  /*
   * Reposition spline and scale dimensionwise
   *
   * Scales a spline along its different dimensions and transposes its position.
   * This can be used as an internal function when a spline is mapped into the
   * unit cube, but is also handy when the same operation is performed on a
   * group of splines where the scaling parameters and transposition vectors
   * stay constant for the entire group.
   *
   * @param transposition PointType describing the first corner of the nd cuboid
   * @param stretch PointType describing the second corner of the nd cuboid
   */
  constexpr BezierSpline& TransposeAndScale(
      const PhysicalPointType& transposition,
      const PhysicalPointType& scale_vector);

  /*
   * Fit into unit cube
   *
   * Takes spline and fits it into the unit cuboid (important for spline
   * composition)
   */
  constexpr BezierSpline& FitToUnitCube();

  /// Friend injection Substraction
  friend constexpr BezierSpline operator-(const PhysicalPointType& point_shift,
                                          const BezierSpline& original_spline) {
    return -(original_spline - point_shift);
  }

  /*
   * Functional Composition between two splines
   *
   * Compose two splines, taking the (*this) spline as the outer funtion and the
   * function argument as the inner function. This works so long as the
   * parametric dimension of the outer function matches the physical dimension
   * of the inner function.
   */
  template <std::size_t parametric_dimension_inner_spline,
            typename PointTypeRHS, typename ScalarRHS>
  constexpr BezierSpline<parametric_dimension_inner_spline, PhysicalPointType,
                         decltype(ScalarType_{} * ScalarRHS{})>
  Compose(const BezierSpline<parametric_dimension_inner_spline, PointTypeRHS,
                             ScalarRHS>& inner_function) const;

  /*
   * Composition between mutliple splines from a spline group
   *
   * Performes a composition between multple splines, which can be used to
   * construct microstructures. After the return group is instantiated, the
   * composition is performed elementwise.
   */
  template <std::size_t parametric_dimension_inner_spline,
            typename PointTypeRHS, typename ScalarRHS>
  constexpr BezierSplineGroup<parametric_dimension_inner_spline,
                              PhysicalPointType,
                              decltype(ScalarType_{} * ScalarRHS{})>
  Compose(
      const BezierSplineGroup<parametric_dimension_inner_spline, PointTypeRHS,
                              ScalarRHS>& inner_function_group) const;

  /*
   * Split the Bezier Spline into two distinct subdivisions
   *
   * Splits the Spline along a specific dimension and returns a group
   * representing the same domain over two splines.
   */
  constexpr BezierSplineGroup<parametric_dimension, PhysicalPointType,
                              ScalarType>
  SplitAtPosition(const ScalarType& splitting_plane,
                  const IndexingType splitting_dimension = 0) const;

  /*
   * Split the Bezier Spline into several subdivisions
   *
   * Splits the Spline along a specific dimension and returns a group
   * representing the same domain over several splines.
   */
  constexpr BezierSplineGroup<parametric_dimension, PhysicalPointType,
                              ScalarType>
  SplitAtPosition(const std::vector<ScalarType>& splitting_planes,
                  const IndexingType splitting_dimension = 0) const;
};

#include "bezierManipulation/src/bezier_spline.inc"

}  // namespace beziermanipulation

#endif  // SRC_BEZIER_SPLINE_HPP
