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

#ifndef SOURCES_VECTORSPACES_VECTOR_SPACE_HPP_
#define SOURCES_VECTORSPACES_VECTOR_SPACE_HPP_

#include <algorithm>
#include <functional>
#include <utility>

#include "Sources/Utilities/error_handling.hpp"
#include "Sources/Utilities/named_type.hpp"
#include "Sources/Utilities/numeric_operations.hpp"
#include "Sources/Utilities/std_container_operations.hpp"
#include "Sources/Utilities/string_operations.hpp"

namespace splinelib::sources::vector_spaces {

template<int dimensionality> class VectorSpace;

template<int dimensionality>
bool IsEqual(VectorSpace<dimensionality> const &lhs, VectorSpace<dimensionality> const &rhs,
             Tolerance const &tolerance = kEpsilon);
template<int dimensionality>
bool operator==(VectorSpace<dimensionality> const &lhs, VectorSpace<dimensionality> const &rhs);

// VectorSpaces group coordinates.
//
// Example:
//   using VectorSpace3d = VectorSpace<3>;
//   using Coordinate = VectorSpace3d::Coordinate_;
//   using ScalarCoordinate = Coordinate::value_type;
//   ScalarCoordinate const k0_0{}, k1_0{1.0};
//   VectorSpace3d const vector_space{{{k0_0, k0_0, k0_0}, {k1_0, k0_0, k0_0}, {k0_0, k1_0, k0_0}, {k1_0, k1_0, k0_0}}};
//   int const &four = vector_space.GetNumberOfCoordinates();
//   Coordinate const &coordinate = vector_space[Index{1}];  // Coordinate P_1 = {1.0, 0.0, 0.0}.
//   ScalarCoordinate const &one_point_zero = vector_space.DetermineMaximumDistanceFromOrigin();
template<int dimensionality>
class VectorSpace {
 public:
  using Coordinate_ = Array<Coordinate, dimensionality>;
  using Coordinates_ = Vector<Coordinate_>;
  using OutputInformation_ = Tuple<Vector<StringArray<dimensionality>>>;

  VectorSpace() = default;
  explicit VectorSpace(Coordinates_ coordinates);
  VectorSpace(VectorSpace const &other) = default;
  VectorSpace(VectorSpace &&other) noexcept = default;
  VectorSpace & operator=(VectorSpace const &rhs) = default;
  VectorSpace & operator=(VectorSpace &&rhs) noexcept = default;
  virtual ~VectorSpace() = default;

  // Comparison based on tolerance.
  friend bool IsEqual<dimensionality>(VectorSpace const &lhs, VectorSpace const &rhs, Tolerance const &tolerance);
  // Comparison based on numeric_operations::GetEpsilon<Tolerance>().
  friend bool operator==<dimensionality>(VectorSpace const &lhs, VectorSpace const &rhs);
  virtual Coordinate_ const & operator[](Index const &coordinate) const;

  virtual int GetNumberOfCoordinates() const;
  virtual void Replace(Index const &coordinate_index, Coordinate_ coordinate);
  virtual void Insert(Index const &coordinate_index, Coordinate_ coordinate);
  virtual void Erase(Index const &coordinate_index);

  virtual Coordinate DetermineMaximumDistanceFromOrigin(Tolerance const &tolerance = kEpsilon) const;
  virtual OutputInformation_ Write(Precision const &precision = kPrecision) const;

 protected:
  Coordinates_ coordinates_;

 private:
#ifndef NDEBUG
  void ThrowIfIndexIsInvalid(Index const &coordinate) const;
#endif
};

#include "Sources/VectorSpaces/vector_space.inc"

}  // namespace splinelib::sources::vector_spaces

#endif  // SOURCES_VECTORSPACES_VECTOR_SPACE_HPP_
