/* Copyright (c) 2018–2021 SplineLib

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit
persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. */

template<int parametric_dimensionality, int dimensionality>
bool IsEqual(Spline<parametric_dimensionality, dimensionality> const &lhs,
             Spline<parametric_dimensionality, dimensionality> const &rhs, Tolerance const &tolerance) {
  using Base = typename Spline<parametric_dimensionality, dimensionality>::Base_;

#ifndef NDEBUG
  try {
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const &exception) {
    Throw(exception, "splinelib::sources::splines::IsEqual::Spline");
  }
#endif
  return ((static_cast<Base const &>(lhs) == static_cast<Base const &>(rhs)) &&
          IsEqual(*lhs.parameter_space_, *rhs.parameter_space_, tolerance));
}

template<int parametric_dimensionality, int dimensionality>
bool operator==(Spline<parametric_dimensionality, dimensionality> const &lhs,
                Spline<parametric_dimensionality, dimensionality> const &rhs) {
  return IsEqual(lhs, rhs);
}

template<int parametric_dimensionality, int dimensionality>
void Spline<parametric_dimensionality, dimensionality>::RefineKnots(Dimension const &dimension, Knots_ knots,
    Multiplicity const &multiplicity, Tolerance const &tolerance) const {
  std::for_each(knots.begin(), knots.end(), [&] (Knot_ const &knot) { InsertKnot(dimension, std::move(knot),
                                                                                 multiplicity, tolerance); });
}

template<int parametric_dimensionality, int dimensionality>
Multiplicity Spline<parametric_dimensionality, dimensionality>::CoarsenKnots(Dimension const &dimension,
    Knots_ const &knots, Tolerance const &tolerance_removal, Multiplicity const &multiplicity,
    Tolerance const &tolerance) const {
  Multiplicity successful_removals{multiplicity};
  std::for_each(knots.begin(), knots.end(), [&] (Knot_ const &knot) {
      successful_removals = std::min(successful_removals, RemoveKnot(dimension, knot, tolerance_removal, multiplicity,
                                                                     tolerance)); });
  return successful_removals;
}

template<int parametric_dimensionality, int dimensionality>
typename Spline<parametric_dimensionality, dimensionality>::Coordinates_
Spline<parametric_dimensionality, dimensionality>::Sample(
    NumberOfParametricCoordinates_ const &number_of_parametric_coordinates, Tolerance const &tolerance) const {
  using ParametricCoordinates = typename ParameterSpace_::ParametricCoordinates_;

  ParametricCoordinates const &parametric_coordinates = parameter_space_->Sample(number_of_parametric_coordinates);
  Coordinates_ coordinates{};
  coordinates.reserve(parametric_coordinates.size());
  std::for_each(parametric_coordinates.begin(), parametric_coordinates.end(),
                [&] (typename ParametricCoordinates::value_type const &parametric_coordinate) {
                    coordinates.emplace_back(operator()(parametric_coordinate, tolerance)); });
  return coordinates;
}

template<int parametric_dimensionality, int dimensionality>
Spline<parametric_dimensionality, dimensionality>::Spline(bool is_rational) :
    SplineItem(parametric_dimensionality, dimensionality, std::move(is_rational)) {}

template<int parametric_dimensionality, int dimensionality>
Spline<parametric_dimensionality, dimensionality>::Spline(SharedPointer<ParameterSpace_> parameter_space,
    bool is_rational) : Spline(std::move(is_rational)) {
  static_assert(parametric_dimensionality > 0, "The parametric dimensionality must be positive");
  static_assert(dimensionality > 0, "The dimensionality must be positive");

  parameter_space_ = std::move(parameter_space);
}

template<int parametric_dimensionality, int dimensionality>
Spline<parametric_dimensionality, dimensionality>::Spline(Spline const &other) : Base_(),
    parameter_space_(std::make_shared<ParameterSpace_>(*other.parameter_space_)) {}

template<int parametric_dimensionality, int dimensionality>
Spline<parametric_dimensionality, dimensionality> &
Spline<parametric_dimensionality, dimensionality>::operator=(Spline const &rhs) {
  Base_::operator=(rhs);
  parameter_space_ = std::make_shared<ParameterSpace_>(*rhs.parameter_space_);
  return *this;
}
