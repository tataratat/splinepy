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

template<int parametric_dimensionality, int geometric_dimensionality, bool is_rational>
auto Generate(NumbersOfNonZeroKnotSpans<parametric_dimensionality> const &numbers_of_non_zero_knot_spans,
              Degrees<parametric_dimensionality> const &degrees) {
  static_assert(parametric_dimensionality <= geometric_dimensionality, "Geometric dimensionality must not be less than "
                                                                       "parametric dimensionality!");

  using DimensionType = Dimension::Type_;
  using IndexType = Index::Type_;
  using ParameterSpace = parameter_spaces::ParameterSpace<parametric_dimensionality>;
  using Coordinates = typename vector_spaces::VectorSpace<geometric_dimensionality>::Coordinates_;
  using Index = typename ParameterSpace::Index_;
  using KnotVectors = typename ParameterSpace::KnotVectors_;
  using std::make_shared;

  KnotVectors knot_vectors;
  for (DimensionType dimension{}; dimension < parametric_dimensionality; ++dimension) {
    using KnotVector = typename KnotVectors::value_type::element_type;
    using Knots = typename KnotVector::Knots_;
    using Knot = typename Knots::value_type;

    IndexType const &number_of_non_zero_knot_spans = numbers_of_non_zero_knot_spans[dimension].Get();
    Knots knots;
    knots.reserve(number_of_non_zero_knot_spans);
    Degree::Type_ const &degree_value = degrees[dimension].Get();
    for (IndexType non_zero_knot_span{}; non_zero_knot_span <= degree_value; ++non_zero_knot_span)
        knots.emplace_back(0.0);
    for (IndexType non_zero_knot_span{0}; non_zero_knot_span < number_of_non_zero_knot_spans; ++non_zero_knot_span)
        knots.emplace_back(static_cast<typename Knot::Type_>(1 + non_zero_knot_span) /
                           number_of_non_zero_knot_spans);
    for (IndexType non_zero_knot_span{}; non_zero_knot_span < degree_value; ++non_zero_knot_span)
        knots.emplace_back(1.0);
    knot_vectors[dimension] = make_shared<KnotVector>(knots);
  }
  SharedPointer<ParameterSpace> parameter_space = make_shared<ParameterSpace>(knot_vectors, degrees);
  typename ParameterSpace::NumberOfBasisFunctions_ const &number_of_basis_functions =
                                                             parameter_space->GetNumberOfBasisFunctions();
  int const &total_number_of_basis_functions = parameter_space->GetTotalNumberOfBasisFunctions();
  Coordinates coordinates;
  coordinates.reserve(total_number_of_basis_functions);
  for (Index index_coordinate{Index::First(number_of_basis_functions)};
       index_coordinate != Index::Behind(number_of_basis_functions); ++index_coordinate) {
    using Coordinate = typename Coordinates::value_type;
    using CoordinateScalar = typename Coordinate::value_type;

    Coordinate coordinate{};
    for (DimensionType dimension{}; dimension < parametric_dimensionality; ++dimension)
        coordinate[dimension] = CoordinateScalar{static_cast<typename CoordinateScalar::Type_>(
            index_coordinate[Dimension{dimension}].Get()) / (number_of_basis_functions[dimension].Get() - 1)};
    coordinates.push_back(coordinate);
  }
  if constexpr (is_rational) {
    using Nurbs = splines::Nurbs<parametric_dimensionality, geometric_dimensionality>;
    using WeightedVectorSpace = typename Nurbs::WeightedVectorSpace_;
    using Weights = typename WeightedVectorSpace::Weights_;

    Weights weights(total_number_of_basis_functions, typename Weights::value_type{2.0});
    return Nurbs{parameter_space, make_shared<WeightedVectorSpace>(coordinates, weights)};
  } else {
    using BSpline = splines::BSpline<parametric_dimensionality, geometric_dimensionality>;

    return BSpline{parameter_space, make_shared<typename BSpline::VectorSpace_>(coordinates)};
  }
}
