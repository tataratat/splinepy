template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr void BezierSpline<parametric_dimension, PhysicalPointType,
                            ScalarType>::UpdateIndexOffsets_() {
  index_offsets[0] = 1;
  for (unsigned int i{1}; i < parametric_dimension; i++)
    index_offsets[i] = index_offsets[i - 1] * (degrees[i - 1] + 1);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr std::array<std::size_t, parametric_dimension> BezierSpline<
    parametric_dimension, PhysicalPointType,
    ScalarType>::LocalToGlobalIndex(const IndexingType& local_index) const {
  std::array<IndexingType, parametric_dimension> indexList{};
  for (unsigned int i{0}; i < parametric_dimension; i++) {
    indexList[i] = (local_index / index_offsets[i]) % (degrees[i] + 1);
  }
  return indexList;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename PhysicalPointLHS, typename ScalarLHS,
          typename PhysicalPointRHS, typename ScalarRHS, typename... T>
constexpr void
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    CombineControlPointsForProduct_(
        const BezierSpline<parametric_dimension, PhysicalPointLHS, ScalarLHS>&
            P_spline,
        const BezierSpline<parametric_dimension, PhysicalPointRHS, ScalarRHS>&
            Q_spline,
        const std::array<IndexingType, parametric_dimension>& ctpsIndex,
        const ScalarType factor, const T&... indices) {
  // Some constant indices and degrees
  const int depth = sizeof...(indices);
  const int k = ctpsIndex[depth];
  const int m = P_spline.GetDegrees()[depth];
  const int n = Q_spline.GetDegrees()[depth];

  // Loop over current parametric domain
  for (int i{std::max(0, k - n)}; i <= std::min(k, m); i++) {
    // Calculate Factor
    const ScalarType lFactor =
        static_cast<ScalarType>(
            utils::FastBinomialCoefficient::choose(m, i) *
            utils::FastBinomialCoefficient::choose(n, k - i)) /
        static_cast<ScalarType>(
            utils::FastBinomialCoefficient::choose(m + n, k));

    // Now decide if continue recursion
    if constexpr ((depth + 1) == parametric_dimension) {
      const std::array<IndexingType, parametric_dimension> ind_lhs{
          static_cast<IndexingType>(indices)..., static_cast<IndexingType>(i)};
      std::array<IndexingType, parametric_dimension> ind_rhs{};
      for (unsigned int j{}; j < parametric_dimension; j++) {
        ind_rhs[j] = ctpsIndex[j] - ind_lhs[j];
      }
      (*this).ControlPoint(ctpsIndex) += P_spline.ControlPoint(ind_lhs) *
                                         Q_spline.ControlPoint(ind_rhs) *
                                         factor * lFactor;
    } else {
      CombineControlPointsForProduct_(P_spline, Q_spline, ctpsIndex,
                                      factor * lFactor, indices..., i);
    }
  }
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename... T>
constexpr const PhysicalPointType&
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::ControlPoint(
    const T... index) const {
  static_assert(sizeof...(T) == parametric_dimension,
                "Unspecified number of indices.");
  unsigned int c_i{0}, i{};
  ((c_i += index_offsets[i++] * index), ...);
  return control_points[c_i];
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename... T>
constexpr PhysicalPointType&
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::ControlPoint(
    const T... index) {
  static_assert(sizeof...(T) == parametric_dimension,
                "Unspecified number of indices.");
  unsigned int c_i{0}, i{};
  ((c_i += index_offsets[i++] * index), ...);
  return control_points[c_i];
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr const PhysicalPointType&
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::ControlPoint(
    const std::array<IndexingType, parametric_dimension>& index) const {
  unsigned int c_i{0}, i{};
  for (unsigned int i{}; i < parametric_dimension; i++) {
    c_i += index_offsets[i] * index[i];
  }
  return control_points[c_i];
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr PhysicalPointType&
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::ControlPoint(
    const std::array<IndexingType, parametric_dimension>& index) {
  unsigned int c_i{0}, i{};
  for (unsigned int i{}; i < parametric_dimension; i++) {
    c_i += index_offsets[i] * index[i];
  }
  return control_points[c_i];
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    OrderElevateAlongParametricDimension(const IndexingType par_dim) {
  // Calculate index Offsets to facilitate working on 1D array
  const unsigned int n_starting_points =
      (NumberOfControlPoints / (degrees[par_dim] + 1));
  const unsigned int starting_point_offset =
      index_offsets[par_dim] * (degrees[par_dim] + 1);
  const int starting_points_per_group = index_offsets[par_dim];
  const int n_groups = n_starting_points / starting_points_per_group;

  // Resize the CTPS vector accordingly
  NumberOfControlPoints =
      NumberOfControlPoints / (degrees[par_dim] + 1) * (degrees[par_dim] + 2);
  control_points.resize(NumberOfControlPoints);
  degrees[par_dim]++;

  // Local Counter
  unsigned int global_index = NumberOfControlPoints - 1;

  // Precalculations
  const ScalarType_ inverse_factor =
      static_cast<ScalarType_>(1) / static_cast<ScalarType_>(degrees[par_dim]);
  const IndexingType variable_offset_factor =
      index_offsets[par_dim] * (degrees[par_dim] - 1);

  // Vector is calculated from back to front, to hinder overwrite
  for (int group_index{n_groups - 1}; group_index >= 0; group_index--) {
    // Local variables
    const unsigned int first_index_in_group =
        group_index * (starting_point_offset);
    IndexingType i = degrees[par_dim] - 1;

    // Fix the last entry for element in the group
    for (int index_in_group{starting_points_per_group - 1}; index_in_group >= 0;
         index_in_group--) {
      control_points[global_index] =
          control_points[first_index_in_group + index_in_group +
                         variable_offset_factor];
      global_index--;
    }

    // Interpolate for all but the first points in the vector, constantly
    // decreasing the counter index algorithm found in
    // https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/Bezier/bezier-elev.html
    // The complex indexation is a result from the row based storage of the
    // control point positions
    for (IndexingType i{degrees[par_dim] - 1}; i > 0; i -= 1) {
      for (int index_in_group{starting_points_per_group - 1};
           index_in_group >= 0; index_in_group--) {
        const ScalarType_ factor = static_cast<ScalarType_>(i) * inverse_factor;

        control_points[global_index] =
            control_points[first_index_in_group + index_in_group +
                           index_offsets[par_dim] * i] *
                (1 - factor) +
            control_points[first_index_in_group + index_in_group +
                           index_offsets[par_dim] * (i - 1)] *
                factor;
        global_index--;
      }
    }

    // Fixate the first entry along each parametric dimension
    for (int index_in_group{starting_points_per_group - 1}; index_in_group >= 0;
         index_in_group--) {
      // No we actually start the algorithm
      control_points[global_index] =
          control_points[first_index_in_group + index_in_group];
      global_index--;
    }
  }
  UpdateIndexOffsets_();
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    DerivativeWRTParametricDimension(const IndexingType par_dim) const {
  // Initialize return value and precompute auxiliary values
  auto new_degrees = degrees;
  new_degrees[par_dim]--;
  BezierSpline derivative{new_degrees};
  const ScalarType factor = static_cast<ScalarType>(degrees[par_dim]);

  // Start looping over all control points in derivative and determine their new
  // values
  for (IndexingType i_point{0}; i_point < derivative.NumberOfControlPoints;
       i_point++) {
    const auto index = derivative.LocalToGlobalIndex(i_point);
    auto next_index = index;
    next_index[par_dim]++;
    derivative.control_points[i_point] =
        factor * (ControlPoint(next_index) - ControlPoint(index));
  }
  return derivative;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr PhysicalPointType
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::Evaluate(
    const PointTypeParametric_& par_coords) const {
  /*
   * Notation:
   * A group the number of all starting points that are subsequent, i.e. the
   * size of a group is the number of subsequent 0s. E.g. in the example:
   * Group size of index 0 is 1, group size of index 1 is 3, etc.. Group
   * offset is thus the distance between two group-starting-points. The idea is
   * to collaps the groups one by one starting from the first parametric
   * dimension going forward. For indexing, see SplitAtPosition(...)
   */

  // Work on copy of ControlPoint
  std::vector<PointTypePhysical_> ctps_copy{control_points};
  IndexingType ctps_to_consider = NumberOfControlPoints;

  for (IndexingType par_dim{0}; par_dim < parametric_dimension; par_dim++) {
    ScalarType_ factor = par_coords[par_dim];
    ScalarType_ inv_factor = 1. - par_coords[par_dim];

    ctps_to_consider /= degrees[par_dim] + 1;

    // For every starting position
    for (IndexingType starting_point_id{0};
         starting_point_id < ctps_to_consider; starting_point_id++) {
      const auto offset =
          (index_offsets[par_dim] * (degrees[par_dim] + 1)) * starting_point_id;
      const auto step_width = index_offsets[par_dim];

      for (IndexingType i{0}; i <= degrees[par_dim]; i++) {
        for (IndexingType j{0}; j < degrees[par_dim] - i; j++) {
          ctps_copy[j * step_width + offset] =
              ctps_copy[j * step_width + offset] * inv_factor +
              ctps_copy[(j + 1) * step_width + offset] * factor;
        }
      }
    }
  }
  return ctps_copy[0];
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr PhysicalPointType BezierSpline<
    parametric_dimension, PhysicalPointType,
    ScalarType>::ForwardEvaluate(const PointTypeParametric_& par_coords) const {
  /*
   * Here we avoid Point operations and try to precompute as many values as
   * possible, this greatly increases the complexity of the necessary code.
   *
   * They are all stored into an array to store them on the stack and avoid
   * allocation
   */
  std::array<std::array<ScalarType, MAX_BINOMIAL_DEGREE + 1>,
             parametric_dimension>
      factors;

  for (IndexingType par_dim{0}; par_dim < parametric_dimension; par_dim++) {
    ScalarType_ cfactor = static_cast<ScalarType>(1.);
    const ScalarType_ t_value = par_coords[par_dim];
    const ScalarType_ inv_t_value = 1. - par_coords[par_dim];
    const IndexingType c_degree = degrees[par_dim];
    // Loop over dimensionwise ctps
    for (IndexingType i_ctps{0}; i_ctps < c_degree + 1; i_ctps++) {
      factors[par_dim][i_ctps] =
          static_cast<ScalarType>(
              utils::FastBinomialCoefficient::choose(c_degree, i_ctps)) *
          cfactor;
      cfactor *= t_value;
    }
    // Backwards assignment
    cfactor = static_cast<ScalarType>(1.);
    for (IndexingType i_ctps{0}; i_ctps < c_degree + 1; i_ctps++) {
      factors[par_dim][c_degree - i_ctps] *= cfactor;
      cfactor *= inv_t_value;
    }
  }

  /*
   * Now that all values have been precomputed we have to sum up the individual
   * contributions using (compile-time) recursion
   */
  PointTypePhysical_ evaluation_point{};

  return AddUpContributionsToControlPointVector_(evaluation_point, factors,
                                                 static_cast<ScalarType>(1.));
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename... Indices>
constexpr PhysicalPointType
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    AddUpContributionsToControlPointVector_(
        PhysicalPointType& evaluation_point,
        const std::array<std::array<ScalarType, MAX_BINOMIAL_DEGREE + 1>,
                         parametric_dimension>& factors,
        const ScalarType& factor_product, const Indices&... indices) const {
  // Check that all are indices

  // Loop over current parametric dimensions
  const IndexingType c_par_dim = sizeof...(Indices);
  const IndexingType i_ctps_per_par_dim = degrees[c_par_dim] + 1;

  for (IndexingType i_ctps{}; i_ctps < i_ctps_per_par_dim; i_ctps++) {
    if constexpr (sizeof...(indices) + 1 == parametric_dimension) {
      evaluation_point += factor_product * factors[c_par_dim][i_ctps] *
                          ControlPoint(indices..., i_ctps);
    } else {
      AddUpContributionsToControlPointVector_(
          evaluation_point, factors,
          factors[c_par_dim][i_ctps] * factor_product, indices..., i_ctps);
    }
  }

  return evaluation_point;
};

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr bool
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator==(
    const BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
        rhs) const {
  // Check if degrees fit and if control_points are the same
  if (rhs.GetDegrees() == degrees) {
    return rhs.control_points == control_points;
  } else {
    return false;
  }
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator*=(
    const ScalarType& scalar) {
  for (std::size_t i{}; i < NumberOfControlPoints; i++) {
    control_points[i] *= scalar;
  }
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator+=(
    const PhysicalPointType& point_shift) {
  for (std::size_t i{}; i < NumberOfControlPoints; i++) {
    control_points[i] += point_shift;
  }
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator+(
    const PhysicalPointType& point_shift) const {
  BezierSpline scaled_spline{(*this)};
  scaled_spline += point_shift;
  return scaled_spline;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator-=(
    const PhysicalPointType& point_shift) {
  for (std::size_t i{}; i < NumberOfControlPoints; i++) {
    control_points[i] -= point_shift;
  }
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator-(
    const PhysicalPointType& point_shift) const {
  BezierSpline scaled_spline{(*this)};
  scaled_spline -= point_shift;
  return scaled_spline;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator-()
    const {
  BezierSpline inverted_spline{(*this)};
  for (IndexingType i{}; i < NumberOfControlPoints; i++) {
    inverted_spline.control_points[i] = -control_points[i];
  }
  return inverted_spline;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator*(
    const ScalarType& scalar) const {
  BezierSpline scaled_spline{(*this)};
  scaled_spline *= scalar;
  return scaled_spline;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename PointTypeRHS, typename ScalarRHS>
constexpr BezierSpline<parametric_dimension,
                       decltype(PhysicalPointType{} + PointTypeRHS{}),
                       decltype(ScalarType{} + ScalarRHS{})>
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator+(
    const BezierSpline<parametric_dimension, PointTypeRHS, ScalarRHS>& rhs)
    const {
  // Initialize return value
  using PointTypeReturnT = decltype(PhysicalPointType{} + PointTypeRHS{});
  using ScalarReturnT = decltype(ScalarType_{} * ScalarRHS{});

  BezierSpline<parametric_dimension, PointTypeReturnT, ScalarReturnT>
      return_spline(degrees, control_points);

  // Check if the right hand side requires a copy as it should not be
  // altered for this purpose
  bool rhs_needs_copy = false;
  for (IndexingType par_dim{0}; par_dim < parametric_dimension; par_dim++) {
    rhs_needs_copy =
        rhs_needs_copy || (degrees[par_dim] > rhs.GetDegrees()[par_dim]);
  }

  // Increase the order of the copied spline to be greater or equal to the
  // RHSs order
  for (IndexingType par_dim{0}; par_dim < parametric_dimension; par_dim++) {
    while (rhs.GetDegrees()[par_dim] > return_spline.GetDegrees()[par_dim]) {
      return_spline.OrderElevateAlongParametricDimension(par_dim);
    }
  }
  if (rhs_needs_copy) {
    // use commutativity of addition to create a copy of rhs
    return rhs + return_spline;
  } else {
    for (IndexingType i_ctps{}; i_ctps < return_spline.NumberOfControlPoints;
         i_ctps++) {
      return_spline.control_points[i_ctps] += rhs.control_points[i_ctps];
    }
    return return_spline;
  }
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator+=(
    BezierSpline<parametric_dimension, PhysicalPointType, ScalarType> rhs) {
  // Increase the order of (*this) spline to be greater or equal to the
  // RHSs order
  for (IndexingType par_dim{0}; par_dim < parametric_dimension; par_dim++) {
    while (rhs.GetDegrees()[par_dim] > degrees[par_dim]) {
      OrderElevateAlongParametricDimension(par_dim);
    }
  }
  // Do the same thing to the rhs
  for (IndexingType par_dim{0}; par_dim < parametric_dimension; par_dim++) {
    while (rhs.GetDegrees()[par_dim] < degrees[par_dim]) {
      rhs.OrderElevateAlongParametricDimension(par_dim);
    }
  }
  // Add control points
  for (IndexingType i_ctps{}; i_ctps < NumberOfControlPoints; i_ctps++) {
    control_points[i_ctps] += rhs.control_points[i_ctps];
  }
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, ScalarType, ScalarType>
BezierSpline<parametric_dimension, PhysicalPointType,
             ScalarType>::ExtractDimension(unsigned int dimension) const {
  if constexpr (std::is_arithmetic_v<PhysicalPointType>) {
    assert(dimension == 0ul);
    return (*this);
  } else {
    assert(dimension < PointTypePhysical_::kSpatialDimension);
    BezierSpline<parametric_dimension, ScalarType, ScalarType> extracted_spline(
        degrees);
    for (std::size_t i{}; i < NumberOfControlPoints; i++) {
      extracted_spline.control_points[i] = control_points[i][dimension];
    }
    return extracted_spline;
  }
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::RaisePower(
    const unsigned int power) const {
  // @TODO use log2(power) algorithm that squares the result to minimize
  // multiplications
  static_assert(std::is_scalar_v<PhysicalPointType>,
                "Only Scalar-type Splines can be raised to a power.");
  assert(("Not implemented, as raising to 0 would be inefficient.", power > 0));
  BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>
      power_spline{(*this)};
  for (int i{1}; i < power; i++) {
    power_spline = power_spline * (*this);
  }
  return power_spline;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr PhysicalPointType
BezierSpline<parametric_dimension, PhysicalPointType,
             ScalarType>::MaximumCorner() const {
  assert(("Spline seems to be non-defined or has no control-points",
          NumberOfControlPoints > 0));
  PhysicalPointType maximum{control_points[0]};
  for (IndexingType i_dim{}; i_dim < PhysicalPointType::kSpatialDimension;
       i_dim++) {
    for (IndexingType i{1}; i < NumberOfControlPoints; i++) {
      maximum[i_dim] = std::max(maximum[i_dim], control_points[i][i_dim]);
    }
  }
  return maximum;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr PhysicalPointType
BezierSpline<parametric_dimension, PhysicalPointType,
             ScalarType>::MinimumCorner() const {
  assert(("Spline seems to be non-defined or has no control-points",
          NumberOfControlPoints > 0));
  PhysicalPointType minimum{control_points[0]};
  for (IndexingType i_dim{}; i_dim < PhysicalPointType::kSpatialDimension;
       i_dim++) {
    for (IndexingType i{1}; i < NumberOfControlPoints; i++) {
      minimum[i_dim] = std::min(minimum[i_dim], control_points[i][i_dim]);
    }
  }
  return minimum;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    TransposeAndScale(const PhysicalPointType& transposition,
                      const PhysicalPointType& scale_vector) {
  for (IndexingType i_dim{}; i_dim < PhysicalPointType::kSpatialDimension;
       i_dim++) {
    for (IndexingType i{}; i < NumberOfControlPoints; i++) {
      control_points[i][i_dim] =
          (transposition[i_dim] + control_points[i][i_dim]) *
          scale_vector[i_dim];
    }
  }
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
BezierSpline<parametric_dimension, PhysicalPointType,
             ScalarType>::FitToUnitCube() {
  const auto maximum_corner = MaximumCorner();
  const auto minimum_corner = MinimumCorner();
  auto scale_vector = maximum_corner - minimum_corner;
  for (IndexingType i_dim{}; i_dim < PhysicalPointType::kSpatialDimension;
       i_dim++) {
    scale_vector[i_dim] = static_cast<ScalarType>(1) / scale_vector[i_dim];
  }
  return TransposeAndScale(-minimum_corner, scale_vector);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr bool BezierSpline<parametric_dimension, PhysicalPointType,
                            ScalarType>::FitsIntoUnitCube() const {
  const auto maximum_corner = MaximumCorner();
  const auto minimum_corner = MinimumCorner();
  for (IndexingType i_dim{}; i_dim < PhysicalPointType::kSpatialDimension;
       i_dim++) {
    if ((maximum_corner[i_dim] > static_cast<ScalarType>(1)) ||
        (minimum_corner[i_dim] < static_cast<ScalarType>(0))) {
      return false;
    }
  }
  return true;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename PointTypeRHS, typename ScalarRHS>
constexpr BezierSpline<parametric_dimension,
                       decltype(PhysicalPointType{} * PointTypeRHS{}),
                       decltype(ScalarType{} * ScalarRHS{})>
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator*(
    const BezierSpline<parametric_dimension, PointTypeRHS,
                       ScalarRHS>& rhs)
    const {  // This multiplication operator is based on the algorithm
  // presented in the thesis from G. Elber (1992)

  // Initialize return value
  using PointTypeReturnT = decltype(PhysicalPointType{} * PointTypeRHS{});
  using ScalarReturnT = decltype(ScalarType{} * ScalarRHS{});

  // Determine the degrees of the resulting spline
  std::array<IndexingType, parametric_dimension> product_degrees;
  for (IndexingType param_dim{}; param_dim < parametric_dimension;
       param_dim++) {
    product_degrees[param_dim] =
        degrees[param_dim] + rhs.GetDegrees()[param_dim];
  }

  // Initialize the return type
  BezierSpline<parametric_dimension, PointTypeReturnT, ScalarReturnT>
      return_spline(product_degrees);

  // Start calculating the new control points
  for (IndexingType i{}; i < return_spline.NumberOfControlPoints; i++) {
    return_spline.CombineControlPointsForProduct_(
        (*this), rhs, return_spline.LocalToGlobalIndex(i),
        static_cast<ScalarReturnT>(1.));
  }
  return return_spline;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <std::size_t parametric_dimension_inner_spline, typename PointTypeRHS,
          typename ScalarRHS>
constexpr BezierSpline<parametric_dimension_inner_spline, PhysicalPointType,
                       decltype(ScalarType{} * ScalarRHS{})>
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::Compose(
    const BezierSpline<parametric_dimension_inner_spline, PointTypeRHS,
                       ScalarRHS>& inner_function) const {
  /// Start the composition with the current spline
  // Initialize return value
  using ScalarReturnT = decltype(ScalarType_{} * ScalarRHS{});
  const IndexingType sum_of_degrees_outer_spline =
      std::accumulate(degrees.begin(), degrees.end(), 0);

  // Check if the inner function fits into the unit cube
  assert(("Inner function must fit unit cube to perform composition.",
          inner_function.FitsIntoUnitCube()));
  // New degrees
  std::array<IndexingType, parametric_dimension_inner_spline> new_degrees{
      inner_function.GetDegrees()};
  for (IndexingType i{}; i < parametric_dimension_inner_spline; i++)
    new_degrees[i] *= sum_of_degrees_outer_spline;

  BezierSpline<parametric_dimension_inner_spline, PhysicalPointType,
               ScalarReturnT>
      composition{new_degrees};

  // Check Dimensions (Scalar Value splines to be treaded as kSpDim = 1)
  if constexpr (std::is_arithmetic_v<PointTypeRHS>) {
    static_assert(parametric_dimension == 1, "Dimension mismatch");
  } else {
    static_assert(PointTypeRHS::kSpatialDimension == parametric_dimension,
                  "Dimension mismatch");
  }

  // Extraction and Precomputation
  /*
   * 1. Extract splines, to get the individual physical coordinates
   *    x_i (xi, eta, ...). These are stored temporarily in a local variable
   * 2. Precompute the "Basis function"-splines, i.e., the spline functions
   *    inserted into the definition of the Bernstein polynomials. This might be
   *    storage inefficient, but reduces the total number of multiplications
   *    required.
   */
  std::array<std::vector<BezierSpline<parametric_dimension_inner_spline,
                                      ScalarReturnT, ScalarReturnT>>,
             parametric_dimension>
      bernstein_basis_splines{};
  for (std::size_t i_outer_parametric_dimension{};
       i_outer_parametric_dimension < parametric_dimension;
       i_outer_parametric_dimension++) {
    // Store current parametric dimensions entry as a reference
    auto& bernstein_spline =
        bernstein_basis_splines[i_outer_parametric_dimension];
    bernstein_spline.resize(degrees[i_outer_parametric_dimension] + 1);

    // Extract the current dimension of the spline and its inverse
    auto extracted_dimension_multiple =
        BezierSpline<parametric_dimension_inner_spline, ScalarReturnT,
                     ScalarReturnT>{
            std::array<std::size_t,
                       parametric_dimension_inner_spline>{},  // initializes
                                                              // to zero
            std::vector<ScalarReturnT>{static_cast<ScalarReturnT>(1.)}};
    auto inverse_extracted_dimension_mulitple = extracted_dimension_multiple;

    const auto extracted_dimension =
        inner_function.ExtractDimension(i_outer_parametric_dimension);

    const auto inverse_extracted_dimension =
        static_cast<ScalarRHS>(1.) - extracted_dimension;

    // Loop over all possible Bernstein basis functions and assign the n-th
    // entry the extracteed dimension to the power n (see "Free Form Surface
    // Analysis Using a Hybrid of Symbolic and Numeric Computation" - Elber -
    // Eq. 6.2)
    // 1. Multiply c(t) ^ i - terms
    bernstein_spline[0] = extracted_dimension_multiple;
    for (std::size_t i_ctps_in_par_dim{1};
         i_ctps_in_par_dim <= degrees[i_outer_parametric_dimension];
         i_ctps_in_par_dim++) {
      extracted_dimension_multiple =
          extracted_dimension_multiple * extracted_dimension;
      bernstein_spline[i_ctps_in_par_dim] = extracted_dimension_multiple;
    }
    // 2. Multiply (1 - c(t)) ^ (n-i) - terms
    for (std::size_t i_ctps_in_par_dim{1};
         i_ctps_in_par_dim <= degrees[i_outer_parametric_dimension];
         i_ctps_in_par_dim++) {
      inverse_extracted_dimension_mulitple =
          inverse_extracted_dimension_mulitple * inverse_extracted_dimension;
      bernstein_spline[degrees[i_outer_parametric_dimension] -
                       i_ctps_in_par_dim] =
          bernstein_spline[degrees[i_outer_parametric_dimension] -
                           i_ctps_in_par_dim] *
          inverse_extracted_dimension_mulitple;
    }
    // 3. Multiply coefficients
    for (std::size_t i_ctps_in_par_dim{1};
         i_ctps_in_par_dim < degrees[i_outer_parametric_dimension];
         i_ctps_in_par_dim++) {
      bernstein_spline[degrees[i_outer_parametric_dimension] -
                       i_ctps_in_par_dim] *=
          utils::FastBinomialCoefficient::choose(
              degrees[i_outer_parametric_dimension], i_ctps_in_par_dim);
    }
  }

  /*
   * @todo The current loop can be made more efficient by precomputing spline
   * products
   *
   * for the three-dimensional case C = T_{ijk} Bi Bj Bk
   * Looping over all i j k (as is currently done) is very inefficient, as many
   * splines Products are ccomputed multiple times (ijk). If the product can be
   * organized as follows, the computational cost can significantly be reduced.
   * for i:
   *    for j:
   *        Bij = Bi * Bj
   *        for k:
   *            Bijk = Bij * Bk
   *
   * In this example for every k -- (k-1) fewer spline products are computed as
   * Bij can be reused for every iteration of k. This effect becomes more
   * important for parametric dimensions of the deformation function
   */
  // Loop over control points of the outer spline to start the interpolation
  for (unsigned int i_control_point_outer_spline{};
       i_control_point_outer_spline < NumberOfControlPoints;
       i_control_point_outer_spline++) {
    // Retrieve indices
    const auto outer_spline_ctps_indices =
        LocalToGlobalIndex(i_control_point_outer_spline);

    // Initialize current factor
    auto bernstein_basis_tensor_product =
        bernstein_basis_splines[0][outer_spline_ctps_indices[0]];
    for (std::size_t i_outer_parametric_dimension{1};
         i_outer_parametric_dimension < parametric_dimension;
         i_outer_parametric_dimension++) {
      bernstein_basis_tensor_product =
          bernstein_basis_tensor_product *
          bernstein_basis_splines
              [i_outer_parametric_dimension]
              [outer_spline_ctps_indices[i_outer_parametric_dimension]];
    }

    // Check if the bernstein_basis_tensor_product has the degrees (and number
    // of control_points) as expected
    assert(("Degree mismatch",
            new_degrees == bernstein_basis_tensor_product.GetDegrees()));

    // Now that all factors have been precomputed we can determine the
    // control-points of the resulting spline
    for (IndexingType i_control_point_composition{};
         i_control_point_composition < composition.NumberOfControlPoints;
         i_control_point_composition++) {
      // Start interpolation
      composition.control_points[i_control_point_composition] +=
          bernstein_basis_tensor_product
              .control_points[i_control_point_composition] *
          control_points[i_control_point_outer_spline];
    }
  }

  return composition;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <std::size_t parametric_dimension_inner_spline, typename PointTypeRHS,
          typename ScalarRHS>
constexpr BezierSplineGroup<parametric_dimension_inner_spline,
                            PhysicalPointType,
                            decltype(ScalarType{} * ScalarRHS{})>
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::Compose(
    const BezierSplineGroup<parametric_dimension_inner_spline, PointTypeRHS,
                            ScalarRHS>& inner_function_group) const {
  /// Start the composition with the current spline
  // Initialize return value
  using ScalarReturnT = decltype(ScalarType_{} * ScalarRHS{});
  const IndexingType number_of_splines{inner_function_group.size()};
  BezierSplineGroup<parametric_dimension_inner_spline, PhysicalPointType,
                    decltype(ScalarType{} * ScalarRHS{})>
      composition_group{number_of_splines};

  // Check if the inner function fits into the unit cube
  assert(("Inner function must fit unit cube to perform composition.",
          inner_function_group.FitsIntoUnitCube()));

  // Check Dimensions
  static_assert(PointTypeRHS::kSpatialDimension == parametric_dimension,
                "Dimension mismatch");

  for (IndexingType i_spline{}; i_spline < number_of_splines; i_spline++) {
    composition_group[i_spline] = Compose(inner_function_group[i_spline]);
  }
  return composition_group;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
BezierSpline<parametric_dimension, PhysicalPointType,
             ScalarType>::SplitAtPosition(const ScalarType& splitting_plane,
                                          const IndexingType
                                              splitting_dimension) const {
  /*
   * Notation:
   * A group the number of all starting points that are subsequent, i.e. the
   * size of a group is the number of subsequent 0s. E.g. in the example:
   * Group size of index 0 is 1, group size of index 1 is 3, etc.. Group
   * offset is thus the distance between two group-starting-points. As before,
   * the index offset is the offset to the next layer with the same
   * corresponding local indices in the other parametric dimensions.
   *
   *    local indices : global index
   *    0   0   0 : 0   -   -
   *    1   0   0 : 1   |   |  group 0 for index 1
   *    2   0   0 : 2   |   -
   *    0   1   0 : 3   |
   *    1   1   0 : 4   | group offset
   *    2   1   0 : 5   |
   *    0   0   1 : 6   -   -
   *    1   0   1 : 7       |  group 1 for index 1
   *    2   0   1 : 8       -
   *    0   1   1 : 9
   *    1   1   1 : 10
   *    2   1   1 : 11
   */

  // Check bounds (as default knot-vector is [0,0,..,1,1,...])
  assert(("Split Plane out of bounds",
          ((splitting_plane > 0.) && (splitting_plane < 1.))));

  // Initialize return values
  BezierSpline first_part{(*this)}, second_part{(*this)};
  BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
      spline_group_split{first_part, second_part};

  // Refernces for simplified access
  auto& ctps_first = spline_group_split[0].control_points;
  auto& ctps_second = spline_group_split[1].control_points;

  // Control points per layer normal to parametric direction
  IndexingType ctps_per_layer =
      NumberOfControlPoints / (degrees[splitting_dimension] + 1);

  // Evaluation factors
  const ScalarType_ factor = splitting_plane;
  const ScalarType_ inv_factor = 1. - splitting_plane;

  // Auxiliary index offsets
  const IndexingType starting_points_per_index_group =
      index_offsets[splitting_dimension];
  const IndexingType group_offset =
      index_offsets[splitting_dimension] * (degrees[splitting_dimension] + 1);
  const IndexingType index_offset = index_offsets[splitting_dimension];

  // Loop over all starting points
  for (IndexingType start_point_id{0}; start_point_id < ctps_per_layer;
       start_point_id++) {
    // Locate position of starting point
    const IndexingType group_id =
        start_point_id / starting_points_per_index_group;
    const IndexingType id_in_group =
        start_point_id % starting_points_per_index_group;
    const IndexingType global_start_id = group_id * group_offset + id_in_group;

    // Going forward on first part and taking last value to second part
    for (IndexingType i{1}; i <= degrees[splitting_dimension]; i++) {
      // Reversed order to prohibit data overwrite
      for (IndexingType j{degrees[splitting_dimension]}; j >= i; j--) {
        ctps_first[global_start_id + j * index_offset] =
            ctps_first[global_start_id + (j - 1) * index_offset] * inv_factor +
            ctps_first[global_start_id + j * index_offset] * factor;
      }
      ctps_second[global_start_id +
                  (degrees[splitting_dimension] - i) * index_offset] =
          ctps_first[global_start_id +
                     degrees[splitting_dimension] * index_offset];
    }
  }
  return spline_group_split;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    SplitAtPosition(const std::vector<ScalarType>& splitting_planes,
                    const IndexingType splitting_dimension) const {
  // Init return value
  BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
      return_value(splitting_planes.size() + 1);

  // Check if the split plane vector has increasing values
  assert(("Split Plane values must be increasing",
          std::is_sorted(splitting_planes.begin(), splitting_planes.end())));

  // If vector is empty return empty group containing only the Spline itself
  if (splitting_planes.size() == 0) {
    return_value[0] = (*this);
    return return_value;
  }

  // Loop over all split planes and assign the first part to the vector,
  // then split the remainng part at the relativ required position
  auto split_spline = SplitAtPosition(splitting_planes[0], splitting_dimension);
  return_value[0] = split_spline[0];
  for (std::size_t i_split{1}; i_split < splitting_planes.size(); i_split++) {
    // Calculate the position of the split position in the remaining spline
    const double split_in_remaining_part =
        (splitting_planes[i_split] - splitting_planes[i_split - 1]) /
        (1. - splitting_planes[i_split - 1]);
    // Split remainder
    split_spline = split_spline[1].SplitAtPosition(split_in_remaining_part,
                                                   splitting_dimension);
    return_value[i_split] = split_spline[0];
  }
  // Assign last section to last element
  return_value[splitting_planes.size()] = split_spline[1];
  return return_value;
}
