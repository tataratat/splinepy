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

#ifndef SOURCES_VECTORSPACES_WEIGHTED_VECTOR_SPACE_HPP_
#define SOURCES_VECTORSPACES_WEIGHTED_VECTOR_SPACE_HPP_

#include <algorithm>
#include <functional>
#include <iterator>
#include <limits>

#include "Sources/Utilities/error_handling.hpp"
#include "Sources/Utilities/named_type.hpp"
#include "Sources/Utilities/numeric_operations.hpp"
#include "Sources/Utilities/std_container_operations.hpp"
#include "Sources/Utilities/string_operations.hpp"
#include "Sources/VectorSpaces/vector_space.hpp"

namespace splinelib::sources::vector_spaces {

template<int dimensionality> class WeightedVectorSpace;

template<int dimensionality>
bool IsEqual(WeightedVectorSpace<dimensionality> const &lhs, WeightedVectorSpace<dimensionality> const &rhs,
             Tolerance const &tolerance = kEpsilon);
template<int dimensionality>
bool operator==(WeightedVectorSpace<dimensionality> const &lhs, WeightedVectorSpace<dimensionality> const &rhs);

// WeightedVectorSpaces store coordinates and weights together using homogeneous, i.e., weighted coordinates.
//
// Example:
//   using WeightedVectorSpace2d = WeightedVectorSpace<2>;
//   Coordinate const k0_0{}, k2_0{2.0};
//   WeightedVectorSpace2d::HomogeneousCoordinate_ const homogeneous_coordinate1{k2_0, k2_0, k2_0};
//   WeightedVectorSpace2d::Coordinate_ const &coordinate1 = WeightedVectorSpace2d::Project(homogeneous_coordinate1);
//   WeightedVectorSpace3d const weighted_vector_space{{{k0_0, k0_0}, coordinate1}, {Weight{1.0}, Weight{2.0}}};
//   WeightedVectorSpace3d::DetermineMaximumDistanceFromOriginAndMinimumWeight_ const &sqrt2_0_and_1_0 =
//       weighted_vector_space.DetermineMaximumDistanceFromOriginAndMinimumWeight();
template<int dimensionality>
class WeightedVectorSpace : public VectorSpace<dimensionality + 1> {
 private:
  using VectorSpace_ = VectorSpace<dimensionality>;

 public:
  using Base_ = VectorSpace<dimensionality + 1>;
  using Coordinate_ = typename VectorSpace_::Coordinate_;
  using Coordinates_ = typename VectorSpace_::Coordinates_;
  using HomogeneousCoordinate_ = typename Base_::Coordinate_;
  using MaximumDistanceFromOriginAndMinimumWeight_ = Tuple<Coordinate, Weight>;
  using OutputInformation_ = Tuple<Vector<StringArray<dimensionality>>, StringVector>;
  using Weights_ = Vector<Weight>;

  WeightedVectorSpace() = default;
  WeightedVectorSpace(Coordinates_ const &coordinates, Weights_ const &weights);
  WeightedVectorSpace(WeightedVectorSpace const &other) = default;
  WeightedVectorSpace(WeightedVectorSpace &&other) noexcept = default;
  WeightedVectorSpace & operator=(WeightedVectorSpace const &rhs) = default;
  WeightedVectorSpace & operator=(WeightedVectorSpace &&rhs) noexcept = default;
  ~WeightedVectorSpace() override = default;

  // Comparison based on tolerance.
  friend bool IsEqual<dimensionality>(WeightedVectorSpace const &lhs, WeightedVectorSpace const &rhs,
                                      Tolerance const &tolerance);
  // Comparison based on numeric_operations::GetEpsilon<Tolerance>().
  friend bool operator==<dimensionality>(WeightedVectorSpace const &lhs, WeightedVectorSpace const &rhs);

  static Coordinate_ Project(HomogeneousCoordinate_ const &homogeneous_coordinate);
  virtual MaximumDistanceFromOriginAndMinimumWeight_ DetermineMaximumDistanceFromOriginAndMinimumWeight(
      Tolerance const &tolerance = kEpsilon) const;

  virtual OutputInformation_ WriteProjected(Precision const &precision = kPrecision) const;

 private:
  using HomogeneousCoordinates_ = typename Base_::Coordinates_;

  HomogeneousCoordinates_ HomogenizeCoordinates(Coordinates_ const &coordinates, Weights_ const &weights) const;
};

#include "Sources/VectorSpaces/weighted_vector_space.inc"

}  // namespace splinelib::sources::vector_spaces

#endif  // SOURCES_VECTORSPACES_WEIGHTED_VECTOR_SPACE_HPP_
