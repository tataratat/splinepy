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

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename... T>
constexpr PhysicalPointType
RationalBezierSpline<parametric_dimension, PhysicalPointType,
                     ScalarType>::ControlPoint(const T... index) const {
  static_assert(sizeof...(T) == parametric_dimension,
                "Unspecified number of indices.");
  const std::array<IndexingType, parametric_dimension> indices{
      (static_cast<IndexingType>(index), ...)};
  return ControlPoint(indices);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr PhysicalPointType
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    ControlPoint(
        const std::array<IndexingType, parametric_dimension>& index) const {
  return weighted_spline_.ControlPoint(index) /
         weight_function_.ControlPoint(index);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename... T>
constexpr PhysicalPointType&
RationalBezierSpline<parametric_dimension, PhysicalPointType,
                     ScalarType>::WeightedControlPoint(const T... index) {
  static_assert(sizeof...(T) == parametric_dimension,
                "Unspecified number of indices.");
  return weighted_spline_.ControlPoint(index...);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename... T>
constexpr const PhysicalPointType&
RationalBezierSpline<parametric_dimension, PhysicalPointType,
                     ScalarType>::WeightedControlPoint(const T... index) const {
  static_assert(sizeof...(T) == parametric_dimension,
                "Unspecified number of indices.");
  return weighted_spline_.ControlPoint(index...);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr PhysicalPointType&
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    WeightedControlPoint(
        const std::array<IndexingType, parametric_dimension>& index) {
  return weighted_spline_.ControlPoint(index);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr const PhysicalPointType&
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    WeightedControlPoint(
        const std::array<IndexingType, parametric_dimension>& index) const {
  return weighted_spline_.ControlPoint(index);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename... T>
constexpr ScalarType&
RationalBezierSpline<parametric_dimension, PhysicalPointType,
                     ScalarType>::Weight(const T... index) {
  static_assert(sizeof...(T) == parametric_dimension,
                "Unspecified number of indices.");
  return weight_function_.ControlPoint(index...);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename... T>
constexpr const ScalarType&
RationalBezierSpline<parametric_dimension, PhysicalPointType,
                     ScalarType>::Weight(const T... index) const {
  static_assert(sizeof...(T) == parametric_dimension,
                "Unspecified number of indices.");
  return weight_function_.ControlPoint(index...);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr ScalarType&
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    Weight(const std::array<IndexingType, parametric_dimension>& index) {
  return weight_function_.ControlPoint(index);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr const ScalarType&
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    Weight(const std::array<IndexingType, parametric_dimension>& index) const {
  return weight_function_.ControlPoint(index);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr PhysicalPointType RationalBezierSpline<
    parametric_dimension, PhysicalPointType,
    ScalarType>::Evaluate(const PointTypeParametric_& par_coords) const {
  return weighted_spline_.Evaluate(par_coords) /
         weight_function_.Evaluate(par_coords);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr std::array<std::vector<ScalarType>, parametric_dimension>
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    BasisFunctions(const PointTypeParametric_& par_coords) const {
  // Evaluiere Basis Functions
  std::array<std::vector<ScalarType>, parametric_dimension> factors =
      weighted_spline_.BasisFunctions(par_coords);

  // Initialize values
  ScalarType weight_eval{};
  weight_function_.AddUpContributionsToControlPointVector_(
      weight_eval, factors, static_cast<ScalarType>(1.));

  // Add up the individual contributions
  const ScalarType inv_weight_function =
      static_cast<ScalarType>(1.) / weight_eval;
  for (IndexingType par_dim{0}; par_dim < parametric_dimension; par_dim++) {
    for (IndexingType i_basis_function{0};
         i_basis_function < factors[par_dim].size(); i_basis_function++) {
      factors[par_dim][i_basis_function] *= inv_weight_function;
    }
  }
  return factors;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr std::array<std::vector<ScalarType>, parametric_dimension>
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    PolynomialBasisFunctions(const PointTypeParametric_& par_coords) const {
  // Pass to member
  return weighted_spline_.BasisFunctions(par_coords);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr PhysicalPointType RationalBezierSpline<
    parametric_dimension, PhysicalPointType,
    ScalarType>::ForwardEvaluate(const PointTypeParametric_& par_coords) const {
  // Evaluiere Basis Functions
  const auto factors = weighted_spline_.BasisFunctions(par_coords);

  // Initialize values
  PhysicalPointType evaluation_point{};
  ScalarType weight_eval{};

  // Add up the individual contributions
  return weighted_spline_.AddUpContributionsToControlPointVector_(
             evaluation_point, factors, static_cast<ScalarType>(1.)) /
         weight_function_.AddUpContributionsToControlPointVector_(
             weight_eval, factors, static_cast<ScalarType>(1.));
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr RationalBezierSpline<parametric_dimension, PhysicalPointType,
                               ScalarType>&
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    OrderElevateAlongParametricDimension(const IndexingType par_dim) {
  weighted_spline_.OrderElevateAlongParametricDimension(par_dim);
  weight_function_.OrderElevateAlongParametricDimension(par_dim);
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr RationalBezierSpline<parametric_dimension, PhysicalPointType,
                               ScalarType>
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    DerivativeWRTParametricDimension(const IndexingType par_dim) const {
  // Implemented using Quotient-rule
  return RationalBezierSpline(
      // u' * w - w' * u
      weighted_spline_.DerivativeWRTParametricDimension(par_dim) *
              weight_function_ -
          weight_function_.DerivativeWRTParametricDimension(par_dim) *
              weighted_spline_,
      // v^2
      weight_function_ * weight_function_);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr RationalBezierSpline<parametric_dimension, PhysicalPointType,
                               ScalarType>&
RationalBezierSpline<parametric_dimension, PhysicalPointType,
                     ScalarType>::operator*=(const ScalarType& scalar) {
  // Multiplication only reflects in the weighted spline
  weighted_spline_ *= scalar;
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr RationalBezierSpline<parametric_dimension, PhysicalPointType,
                               ScalarType>
RationalBezierSpline<parametric_dimension, PhysicalPointType,
                     ScalarType>::operator*(const ScalarType& scalar) const {
  // Multiplication only reflects in the weighted spline
  return RationalBezierSpline(weighted_spline_ * scalar, weight_function_);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr RationalBezierSpline<parametric_dimension, PhysicalPointType,
                               ScalarType>&
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
operator+=(const PhysicalPointType& point_shift) {
  // Add up contributions to control points with weights
  for (IndexingType i_ctps{}; i_ctps < weighted_spline_.NumberOfControlPoints;
       i_ctps++) {
    weighted_spline_.control_points[i_ctps] +=
        point_shift * weight_function_.control_points[i_ctps];
  }
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr RationalBezierSpline<parametric_dimension, PhysicalPointType,
                               ScalarType>
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
operator+(const PhysicalPointType& point_shift) const {
  // Pass to operator +=
  RationalBezierSpline return_spline{weighted_spline_, weight_function_};
  return_spline += point_shift;
  return return_spline;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr RationalBezierSpline<parametric_dimension, PhysicalPointType,
                               ScalarType>&
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
operator-=(const PhysicalPointType& point_shift) {
  // Add up contributions to control points with weights
  for (IndexingType i_ctps{}; i_ctps < weighted_spline_.NumberOfControlPoints;
       i_ctps++) {
    weighted_spline_.control_points[i_ctps] -=
        point_shift * weight_function_.control_points[i_ctps];
  }
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr RationalBezierSpline<parametric_dimension, PhysicalPointType,
                               ScalarType>
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
operator-(const PhysicalPointType& point_shift) const {
  // Pass to operator +=
  RationalBezierSpline return_spline{weighted_spline_, weight_function_};
  return_spline -= point_shift;
  return return_spline;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr RationalBezierSpline<parametric_dimension, PhysicalPointType,
                               ScalarType>
RationalBezierSpline<parametric_dimension, PhysicalPointType,
                     ScalarType>::operator-() const {
  // Pass to operator +=
  RationalBezierSpline return_spline{-weighted_spline_, weight_function_};
  return return_spline;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename PointTypeRHS, typename ScalarRHS>
constexpr RationalBezierSpline<parametric_dimension,
                               decltype(PhysicalPointType{} + PointTypeRHS{}),
                               decltype(ScalarType{} * ScalarRHS{})>
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
operator+(const RationalBezierSpline<parametric_dimension, PointTypeRHS,
                                     ScalarRHS>& rhs) const {
  if (weight_function_ == rhs.weight_function_) {
    return RationalBezierSpline{weighted_spline_ + rhs.weighted_spline_,
                                weight_function_};
  } else {
    return RationalBezierSpline{weighted_spline_ * rhs.weight_function_ +
                                    weight_function_ * rhs.weighted_spline_,
                                weight_function_ * rhs.weight_function_};
  }
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr RationalBezierSpline<parametric_dimension, PhysicalPointType,
                               ScalarType>&
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
operator+=(const RationalBezierSpline<parametric_dimension, PhysicalPointType,
                                      ScalarType>& rhs) {
  if (weight_function_ == rhs.weight_function_) {
    weighted_spline_ += rhs.weighted_spline_;
    return (*this);
  } else {
    weighted_spline_ = weighted_spline_ * rhs.weight_function_;
    weighted_spline_ += weight_function_ * rhs.weighted_spline_;
    weight_function_ = weight_function_ * rhs.weight_function_;
    return (*this);
  }
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename PointTypeRHS, typename ScalarRHS>
constexpr RationalBezierSpline<parametric_dimension,
                               decltype(PhysicalPointType{} - PointTypeRHS{}),
                               decltype(ScalarType{} * ScalarRHS{})>
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
operator-(const RationalBezierSpline<parametric_dimension, PointTypeRHS,
                                     ScalarRHS>& rhs) const {
  if (weight_function_ == rhs.weight_function_) {
    return RationalBezierSpline{weighted_spline_ - rhs.weighted_spline_,
                                weight_function_};
  } else {
    return RationalBezierSpline{weighted_spline_ * rhs.weight_function_ -
                                    weight_function_ * rhs.weighted_spline_,
                                weight_function_ * rhs.weight_function_};
  }
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr RationalBezierSpline<parametric_dimension, PhysicalPointType,
                               ScalarType>&
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
operator-=(const RationalBezierSpline<parametric_dimension, PhysicalPointType,
                                      ScalarType>& rhs) {
  if (weight_function_ == rhs.weight_function_) {
    weighted_spline_ -= rhs.weighted_spline_;
    return (*this);
  } else {
    weighted_spline_ = weighted_spline_ * rhs.weight_function_;
    weighted_spline_ -= weight_function_ * rhs.weighted_spline_;
    weight_function_ = weight_function_ * rhs.weight_function_;
    return (*this);
  }
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename PointTypeRHS, typename ScalarRHS>
constexpr RationalBezierSpline<parametric_dimension,
                               decltype(PhysicalPointType{} * PointTypeRHS{}),
                               decltype(ScalarType{} * ScalarRHS{})>
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
operator*(const RationalBezierSpline<parametric_dimension, PointTypeRHS,
                                     ScalarRHS>& rhs) const {
  return RationalBezierSpline<parametric_dimension,
                              decltype(PhysicalPointType{} * PointTypeRHS{}),
                              decltype(ScalarType{} * ScalarRHS{})>{
      weighted_spline_ * rhs.weighted_spline_,
      weight_function_ * rhs.weight_function_};
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, ScalarType, ScalarType>
RationalBezierSpline<parametric_dimension, PhysicalPointType,
                     ScalarType>::ExtractDimension(const std::size_t& dimension)
    const {
  return weighted_spline_.ExtractDimension(dimension);
}
