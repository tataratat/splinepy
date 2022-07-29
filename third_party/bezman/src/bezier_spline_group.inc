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

#ifndef SRC_BEZIER_SPLINE_GROUP_INC
#define SRC_BEZIER_SPLINE_GROUP_INC

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr bool BezierSplineGroup<parametric_dimension, PhysicalPointType,
                                 ScalarType>::FitsIntoUnitCube() const {
  for (auto spline : (*this)) {
    if (!spline.FitsIntoUnitCube()) {
      utils::Logger::Logging("Group doesn't fit into the unit cube");
      return false;
    }
  }
  utils::Logger::Logging("Group fits into the unit cube");
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
  utils::Logger::ExtendedInformation("Found max. corner at position:" +
                                     maximum_corner.toString());
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
  utils::Logger::ExtendedInformation("Found max. corner at position:" +
                                     minimum_corner.toString());

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

  utils::Logger::ExtendedInformation(
      "Composing " + std::to_string(number_of_splines_in_group) +
      " splines into the outer function spline");

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
  utils::Logger::Logging("Iterating over spline group for derivatives");
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
                  ScalarType>::ExtractDimension(unsigned int physical_dimension)
    const {
  utils::Logger::ExtendedInformation(
      "Extract component (dimension: " + std::to_string(physical_dimension) +
      ") from spline group");
  // Initialize return value with size of current group
  BezierSplineGroup<parametric_dimension, ScalarType, ScalarType> derived_group{
      this->size()};
  for (IndexingType i_spline{}; i_spline < this->size(); i_spline++) {
    derived_group[i_spline] =
        (*this)[i_spline].ExtractDimension(physical_dimension);
  }
  return derived_group;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
    &BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>::
        AddComponentwise(const BezierSplineGroup &rhs) {
  utils::Logger::Logging(
      "Adding two Bezier Spline Groups comonentwise "
      "to the current Group");

  if (this->size() != rhs.size()) {
    utils::Logger::TerminatingError(
        "Size Mismatch for entry-wise Addition between Spline Groups");
  }
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
  utils::Logger::Logging(
      "Adding two Bezier Spline Groups comonentwise "
      "to the current Group");

  if (this->size() != rhs.size()) {
    utils::Logger::TerminatingError(
        "Size Mismatch for Entry-wise Multiplication between two groups");
  }
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
  (*this).insert(this->end(), rhs.begin(), rhs.end());
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
