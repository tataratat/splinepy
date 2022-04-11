#ifndef SRC_BEZIER_SPLINE_GROUP_INC
#define SRC_BEZIER_SPLINE_GROUP_INC

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr bool BezierSplineGroup<parametric_dimension, PhysicalPointType,
                                 ScalarType>::FitsIntoUnitCube() const {
  for (auto spline : (*this)) {
    if (!spline.FitsIntoUnitCube()) {
      return false;
    }
  }
  return true;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr PhysicalPointType
BezierSplineGroup<parametric_dimension, PhysicalPointType,
                  ScalarType>::MaximumCorner() const {
  PhysicalPointType maximum_corner = (*this)[0].MaximumCorner();
  for (IndexingType i{1}; i < this->size(); i++) {
    const auto maximum_corner_i = (*this)[i].MaximumCorner();
    for (IndexingType i_dim{}; i_dim < PhysicalPointType::kSpatialDimension;
         i_dim++) {
      maximum_corner[i_dim] =
          std::max(maximum_corner[i_dim], maximum_corner_i[i_dim]);
    }
  }
  return maximum_corner;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr PhysicalPointType
BezierSplineGroup<parametric_dimension, PhysicalPointType,
                  ScalarType>::MinimumCorner() const {
  PhysicalPointType minimum_corner = (*this)[0].MinimumCorner();
  for (IndexingType i{1}; i < this->size(); i++) {
    const auto minimum_corner_i = (*this)[i].MinimumCorner();
    for (IndexingType i_dim{}; i_dim < PhysicalPointType::kSpatialDimension;
         i_dim++) {
      minimum_corner[i_dim] =
          std::min(minimum_corner[i_dim], minimum_corner_i[i_dim]);
    }
  }
  return minimum_corner;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
    &BezierSplineGroup<parametric_dimension, PhysicalPointType,
                       ScalarType>::FitToUnitCube() {
  const auto maximum_corner = MaximumCorner();
  const auto minimum_corner = MinimumCorner();
  auto scale_vector = maximum_corner - minimum_corner;
  for (IndexingType i_dim{}; i_dim < PhysicalPointType::kSpatialDimension;
       i_dim++) {
    scale_vector[i_dim] = static_cast<ScalarType>(1) / scale_vector[i_dim];
  }
  for (IndexingType i{1}; i < this->size(); i++) {
    (*this)[i].TransposeAndScale(-minimum_corner, scale_vector);
  }
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>::Compose(
    const SplineBaseType &inner_function) const {
  // Initiate return value
  const auto number_of_splines_in_group = this->size();
  BezierSplineGroup composed_group{number_of_splines_in_group};

  // Compose each part of the structure individually
  for (IndexingType i{}; i < number_of_splines_in_group; i++) {
    composed_group[i] = (*this)[i].Compose(inner_function);
  }
  return composed_group;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>::Compose(
    const BezierSplineGroup &inner_function_group) const {
  // Initiate return value
  const auto number_of_splines_in_outer_group = this->size();
  const auto number_of_splines_in_inner_group = inner_function_group.size();

  BezierSplineGroup composed_group{number_of_splines_in_outer_group *
                                   number_of_splines_in_inner_group};

  // Compose each part of the structure individually
  for (IndexingType i{}; i < number_of_splines_in_inner_group; i++) {
    const auto i_composition = Compose(inner_function_group[i]);
    for (IndexingType j{}; j < number_of_splines_in_outer_group; j++) {
      composed_group[j * number_of_splines_in_inner_group + i] =
          i_composition[j];
    }
  }
  return composed_group;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>::
    DerivativeWRTParametricDimension(const IndexingType par_dim) const {
  // Initialize return value with size of current group
  BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
      derived_group{this->size()};
  for (IndexingType i_spline{}; i_spline < this->size(); i_spline++) {
    derived_group[i_spline] =
        (*this)[i_spline].DerivativeWRTParametricDimension(par_dim);
  }
  return derived_group;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSplineGroup<parametric_dimension, ScalarType, ScalarType>
BezierSplineGroup<parametric_dimension, PhysicalPointType,
                  ScalarType>::ExtractDimension(unsigned int par_dim) const {
  // Initialize return value with size of current group
  BezierSplineGroup<parametric_dimension, ScalarType, ScalarType> derived_group{
      this->size()};
  for (IndexingType i_spline{}; i_spline < this->size(); i_spline++) {
    derived_group[i_spline] = (*this)[i_spline].ExtractDimension(par_dim);
  }
  return derived_group;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
    &BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>::
        AddComponentwise(const BezierSplineGroup &rhs) {
  assert(("Size Mismatch", this->size() == rhs.size()));
  for (IndexingType i_spline{}; i_spline < this->size(); i_spline++) {
    (*this)[i_spline] += rhs[i_spline];
  }
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename PointTypeRHS, typename ScalarRHS>
constexpr BezierSplineGroup<parametric_dimension,
                            decltype(PhysicalPointType{} * PointTypeRHS{}),
                            decltype(ScalarType{} * ScalarRHS{})>
BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>::
    MultiplyComponentwise(
        const BezierSplineGroup<parametric_dimension, PointTypeRHS, ScalarRHS>
            &rhs) const {
  assert(("Size Mismatch", this->size() == rhs.size()));
  BezierSplineGroup<parametric_dimension,
                    decltype(PhysicalPointType{} * PointTypeRHS{}),
                    decltype(ScalarType{} * ScalarRHS{})>
      return_group{this->size()};
  for (IndexingType i_spline{}; i_spline < this->size(); i_spline++) {
    return_group[i_spline] = (*this)[i_spline] * rhs[i_spline];
  }
  return return_group;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
BezierSplineGroup<parametric_dimension, PhysicalPointType,
                  ScalarType>::operator+(const BezierSplineGroup &rhs) const {
  BezierSplineGroup combination{(*this)};
  combination += rhs;
  return combination;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
    &BezierSplineGroup<parametric_dimension, PhysicalPointType,
                       ScalarType>::operator+=(const BezierSplineGroup &rhs) {
  (*this).insert(this->begin(), rhs.begin(), rhs.end());
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
BezierSplineGroup<parametric_dimension, PhysicalPointType,
                  ScalarType>::operator+(const PhysicalPointType &translation)
    const {
  BezierSplineGroup combination{(*this)};
  combination += translation;
  return combination;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
    &BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>::
    operator+=(const PhysicalPointType &translation) {
  for (IndexingType i_spline{}; i_spline < this->size(); i_spline++) {
    (*this)[i_spline] += translation;
  }
  return (*this);
}

#endif  // SRC_BEZIER_SPLINE_GROUP_INC
