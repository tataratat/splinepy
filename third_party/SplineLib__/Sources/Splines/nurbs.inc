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
Nurbs<parametric_dimensionality, dimensionality>::Nurbs() : Base_(true) {}

template<int parametric_dimensionality, int dimensionality>
Nurbs<parametric_dimensionality, dimensionality>::Nurbs(SharedPointer<ParameterSpace_> parameter_space,
    SharedPointer<WeightedVectorSpace_> weighted_vector_space) : Base_(std::move(parameter_space), true) {
  using std::to_string;

#ifndef NDEBUG
  int const &total_number_of_basis_functions = Base_::parameter_space_->GetTotalNumberOfBasisFunctions(),
            &number_of_coordinates = weighted_vector_space->GetNumberOfCoordinates();
  if (number_of_coordinates != total_number_of_basis_functions)
      Throw(DomainError(to_string(number_of_coordinates) + " coordinates were provided but " +
                        to_string(total_number_of_basis_functions) + " are needed to associate each basis function "
                        "with a weighted coordinate."), "splinelib::src::spl::Nurbs::Nurbs");
#endif
  homogeneous_b_spline_ = std::make_shared<HomogeneousBSpline_>(Base_::parameter_space_, weighted_vector_space);
  weighted_vector_space_ = std::move(weighted_vector_space);
}

template<int parametric_dimensionality, int dimensionality>
Nurbs<parametric_dimensionality, dimensionality>::Nurbs(Nurbs const &other) : Base_(other),
    homogeneous_b_spline_{std::make_shared<HomogeneousBSpline_>(*other.homogeneous_b_spline_)} {}

template<int parametric_dimensionality, int dimensionality>
Nurbs<parametric_dimensionality, dimensionality> &
Nurbs<parametric_dimensionality, dimensionality>::operator=(Nurbs const &rhs) {
  Base_::operator=(rhs);
  homogeneous_b_spline_ = std::make_shared<HomogeneousBSpline_>(*rhs.homogeneous_b_spline_);
  return *this;
}

template<int parametric_dimensionality, int dimensionality>
bool IsEqual(Nurbs<parametric_dimensionality, dimensionality> const &lhs,
             Nurbs<parametric_dimensionality, dimensionality> const &rhs, Tolerance const &tolerance) {
  return IsEqual(*lhs.homogeneous_b_spline_, *rhs.homogeneous_b_spline_, tolerance);
}

template<int parametric_dimensionality, int dimensionality>
bool operator==(Nurbs<parametric_dimensionality, dimensionality> const &lhs,
                Nurbs<parametric_dimensionality, dimensionality> const &rhs) {
  return IsEqual(lhs, rhs);
}

template<int parametric_dimensionality, int dimensionality>
typename Spline<parametric_dimensionality, dimensionality>::Coordinate_
Nurbs<parametric_dimensionality, dimensionality>::operator()(ParametricCoordinate_ const &parametric_coordinate,
                                                             Tolerance const &tolerance) const {
#ifndef NDEBUG
  try {
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const &exception) { Throw(exception, "splinelib::sources::splines::Nurbs::operator()"); }
#endif
  return WeightedVectorSpace_::Project((*homogeneous_b_spline_)(parametric_coordinate));
}

// See NURBS book Eq. (4.20) (cf. Generalized Leibniz rule at <https://en.wikipedia.org/wiki/General_Leibniz_rule>).
template<int parametric_dimensionality, int dimensionality>
typename Spline<parametric_dimensionality, dimensionality>::Coordinate_
Nurbs<parametric_dimensionality, dimensionality>::operator()(ParametricCoordinate_ const &parametric_coordinate,
    Derivative_ const &derivative, Tolerance const &tolerance) const {
  using BinomialCoefficients = Vector<Derivative::Type_>;
  using DemandForPartialDerivatives = std::deque<bool>;  // std::vector<bool> is not a proper STD container.
  using Index = typename Base_::Index_;
  using IndexLength = typename Index::Length_;
  using IndexValue = typename Index::Value_;
  using ScalarIndexValue = typename IndexValue::value_type;
  using ScalarIndexValueType = typename ScalarIndexValue::Type_;
  using utilities::std_container_operations::TransformNamedTypes;

#ifndef NDEBUG
  try {
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const &exception) { Throw(exception, "splinelib::sources::splines::Nurbs::operator()"); }
#endif
  IndexLength derivative_length;
  std::transform(derivative.begin(), derivative.end(), derivative_length.begin(),
      [] (Derivative const &derivative_for_dimension) { return Length{derivative_for_dimension.Get() + 1}; });
  Index const &zero_derivative = Index::First(derivative_length), &end_derivative = Index::Behind(derivative_length);
  int const &total_number_of_derivatives = zero_derivative.GetTotalNumberOfIndices();

  // Evaluate all homogeneous derivatives of same or lower order.
  typename HomogeneousBSpline_::VectorSpace_::Coordinates_ homogeneous_b_spline_derivatives;
  homogeneous_b_spline_derivatives.reserve(total_number_of_derivatives);
  for (Index current_derivative{zero_derivative}; current_derivative != end_derivative; ++current_derivative)
      homogeneous_b_spline_derivatives.push_back((*homogeneous_b_spline_)(parametric_coordinate,
                                                     TransformNamedTypes<Derivative_>(current_derivative.GetIndex())));

  // Compute binomial coefficients of derivatives up to maximum order.
  Vector<BinomialCoefficients> binomial_coefficients;
  ScalarIndexValueType const &number_of_derivatives = ((*std::max_element(derivative.begin(), derivative.end())).Get() +
                                                       1);
  binomial_coefficients.reserve(number_of_derivatives);
  Derivative::ForEach(0, number_of_derivatives, [&] (Derivative const &current_derivative) {
      Derivative::Type_ const &current_derivative_value = current_derivative.Get(),
                              &current_number_of_derivatives = (current_derivative_value + 1);
      BinomialCoefficients current_binomial_coefficients;
      current_binomial_coefficients.reserve(current_number_of_derivatives);
      Derivative::ForEach(0, current_number_of_derivatives, [&] (Derivative const &lower_derivative) {
          current_binomial_coefficients.push_back(utilities::math_operations::ComputeBinomialCoefficient(
                                                      current_derivative_value, lower_derivative.Get())); });
      binomial_coefficients.push_back(current_binomial_coefficients); });

  // Evaluate requested rational derivative requiring evaluation of all rational derivatives of lower order.
  typename WeightedVectorSpace_::Coordinates_ rational_derivatives;
  rational_derivatives.reserve(total_number_of_derivatives);
  for (Index current_derivative{zero_derivative}; current_derivative != end_derivative; ++current_derivative) {
    IndexValue const &current_derivative_value = current_derivative.GetIndex();
    typename HomogeneousBSpline_::Coordinate_ const &homogeneous_b_spline_derivative =
        homogeneous_b_spline_derivatives[current_derivative.GetIndex1d().Get()];
    Coordinate_ homogeneous_derivative;
    std::copy(homogeneous_b_spline_derivative.begin(), std::prev(homogeneous_b_spline_derivative.end()),
              homogeneous_derivative.begin());
    // Subtract contributions of lower order rational derivatives from current homogeneous derivative by subtracting
    // all incorporated sums.
    DemandForPartialDerivatives are_required_derivatives_currently_required(std::count_if(
        current_derivative_value.begin(), current_derivative_value.end(),
        [] (ScalarIndexValue const &derivative) { return (derivative > ScalarIndexValue{} ? true : false); }), true);
    DemandForPartialDerivatives::iterator const &are_required_derivatives_currently_required_begin =
        are_required_derivatives_currently_required.begin(), &are_required_derivatives_currently_required_end =
            are_required_derivatives_currently_required.end();
    Dimension::ForEach(0, parametric_dimensionality, [&] (Dimension const &dimension) {
        // If all permutations corresponding to the current number of parametric dimensions participating in the sums
        // have been considered, reduce this number by one until all sums have been considered.
        if (std::find(are_required_derivatives_currently_required_begin,
                      are_required_derivatives_currently_required_end, true) !=
                          are_required_derivatives_currently_required_end) {
          do {
            IndexLength lower_derivative_length{utilities::std_container_operations::TransformNamedTypes<IndexLength>(
                                                    current_derivative_value)};
            ScalarIndexValueType permutation_element{};
            Dimension::ForEach(0, parametric_dimensionality, [&] (Dimension const &current_dimension) {
                Dimension::Type_ const &current_dimension_value = current_dimension.Get();
                if (lower_derivative_length[current_dimension_value] > Length{})
                    if (!are_required_derivatives_currently_required[permutation_element++])
                        lower_derivative_length[current_dimension_value] = Length{}; });
            for (Index lower_derivative = Index::First(lower_derivative_length);
                 lower_derivative != Index::Behind(lower_derivative_length); ++lower_derivative) {
              Index const current_lower_derivative{derivative_length, utilities::std_container_operations::Add(
                  utilities::std_container_operations::Subtract(current_derivative_value,
                      TransformNamedTypes<IndexValue>(lower_derivative_length)), lower_derivative.GetIndex())};
              Index const &current_complementary_derivative = (current_derivative -
                                                               current_lower_derivative.GetIndex());
              Degree::Type_ product_of_binomial_coefficients{1};
              Dimension::ForEach(0, parametric_dimensionality, [&] (Dimension const &current_dimension) {
                  product_of_binomial_coefficients *=
                      binomial_coefficients[lower_derivative_length[current_dimension.Get()].Get()][
                                               current_complementary_derivative[current_dimension].Get()]; });
              utilities::std_container_operations::SubtractAndAssignToFirst(homogeneous_derivative,
                  utilities::std_container_operations::Multiply(rational_derivatives[
                      current_lower_derivative.GetIndex1d().Get()], product_of_binomial_coefficients *
                          homogeneous_b_spline_derivatives[current_complementary_derivative.GetIndex1d().Get()]
                                                          [dimensionality].Get()));
            }
          } while (std::next_permutation(are_required_derivatives_currently_required_begin,
                                         are_required_derivatives_currently_required_end));
          are_required_derivatives_currently_required[dimension.Get()] = false;
        } });
    rational_derivatives.push_back(utilities::std_container_operations::Divide(homogeneous_derivative,
                                       homogeneous_b_spline_derivatives[0][dimensionality].Get()));
  }
  return rational_derivatives[total_number_of_derivatives - 1];
}

template<int parametric_dimensionality, int dimensionality>
void Nurbs<parametric_dimensionality, dimensionality>::InsertKnot(Dimension const &dimension, Knot_ knot,
    Multiplicity const &multiplicity, Tolerance const &tolerance) const {
#ifndef NDEBUG
  Message const kName{"splinelib::sources::splines::Nurbs::InsertKnot"};

  try {
    Dimension::ThrowIfNamedIntegerIsOutOfBounds(dimension, parametric_dimensionality - 1);
    Base_::parameter_space_->ThrowIfParametricCoordinateIsOutsideScope(dimension, knot, tolerance);
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (DomainError const &exception) { Throw(exception, kName); }
    catch (OutOfRange const &exception) { Throw(exception, kName); }
    catch (InvalidArgument const &exception) { Throw(exception, kName); }
#endif
  homogeneous_b_spline_->InsertKnot(dimension, knot, multiplicity, tolerance);
}

template<int parametric_dimensionality, int dimensionality>
Multiplicity Nurbs<parametric_dimensionality, dimensionality>::RemoveKnot(Dimension const &dimension,
    Knot_ const &knot, Tolerance const &tolerance_removal, Multiplicity const &multiplicity, Tolerance const &tolerance)
    const {
  using std::get;

#ifndef NDEBUG
  using utilities::numeric_operations::ThrowIfToleranceIsNegative;

  Message const kName{"splinelib::sources::splines::Nurbs::RemoveKnot"};

  try {
    Dimension::ThrowIfNamedIntegerIsOutOfBounds(dimension, parametric_dimensionality - 1);
    ThrowIfToleranceIsNegative(tolerance_removal);
    ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const &exception) { Throw(exception, kName); }
    catch (OutOfRange const &exception) { Throw(exception, kName); }
#endif
  auto const &[maximum_distance_from_origin, minimum_weight] =
      weighted_vector_space_->DetermineMaximumDistanceFromOriginAndMinimumWeight();
  // See NURBS book Eq. (5.30).
  return homogeneous_b_spline_->RemoveKnot(dimension, knot,
      tolerance_removal * (minimum_weight.Get() / (1.0 + maximum_distance_from_origin.Get())), multiplicity, tolerance);
}

template<int parametric_dimensionality, int dimensionality>
void Nurbs<parametric_dimensionality, dimensionality>::ElevateDegree(Dimension const &dimension,
    Multiplicity const &multiplicity, Tolerance const &tolerance) const {
#ifndef NDEBUG
  Message const kName{"splinelib::sources::splines::Nurbs::ElevateDegree"};

  try {
    Dimension::ThrowIfNamedIntegerIsOutOfBounds(dimension, parametric_dimensionality - 1);
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const &exception) { Throw(exception, kName); }
    catch (OutOfRange const &exception) { Throw(exception, kName); }
#endif
  homogeneous_b_spline_->ElevateDegree(dimension, multiplicity, tolerance);
}

template<int parametric_dimensionality, int dimensionality>
bool Nurbs<parametric_dimensionality, dimensionality>::ReduceDegree(Dimension const &dimension,
    Tolerance const &tolerance_removal, Multiplicity const &multiplicity, Tolerance const &tolerance) const {
#ifndef NDEBUG
  using utilities::numeric_operations::ThrowIfToleranceIsNegative;

  Message const kName{"splinelib::sources::splines::Nurbs::ReduceDegree"};

  try {
    Dimension::ThrowIfNamedIntegerIsOutOfBounds(dimension, parametric_dimensionality - 1);
    ThrowIfToleranceIsNegative(tolerance_removal);
    ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const &exception) { Throw(exception, kName); }
    catch (OutOfRange const &exception) { Throw(exception, kName); }
#endif
  return homogeneous_b_spline_->ReduceDegree(dimension, tolerance_removal, multiplicity, tolerance);
}

template<int parametric_dimensionality, int dimensionality>
Coordinate Nurbs<parametric_dimensionality, dimensionality>::ComputeUpperBoundForMaximumDistanceFromOrigin(
    Tolerance const &tolerance) const {
  return std::get<0>(weighted_vector_space_->DetermineMaximumDistanceFromOriginAndMinimumWeight(tolerance));
}

template<int parametric_dimensionality, int dimensionality>
typename Nurbs<parametric_dimensionality, dimensionality>::OutputInformation_
Nurbs<parametric_dimensionality, dimensionality>::Write(Precision const &precision) const {
  return OutputInformation_{Base_::parameter_space_->Write(precision),
                            weighted_vector_space_->WriteProjected(precision)};
}
