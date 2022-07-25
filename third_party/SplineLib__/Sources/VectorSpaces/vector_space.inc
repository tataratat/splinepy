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

template<int dimensionality>
VectorSpace<dimensionality>::VectorSpace(Coordinates_ coordinates) : coordinates_(std::move(coordinates)) {}

template<int dimensionality>
bool IsEqual(VectorSpace<dimensionality> const &lhs, VectorSpace<dimensionality> const &rhs, Tolerance const &tolerance)
{
#ifndef NDEBUG
  try {
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const &exception) {
    Throw(exception, "splinelib::sources::vector_spaces::IsEqual::VectorSpace");
  }
#endif
  return std::equal(lhs.coordinates_.begin(), lhs.coordinates_.end(), rhs.coordinates_.begin(), rhs.coordinates_.end(),
      std::bind(utilities::std_container_operations::DoesContainEqualValues<
          typename VectorSpace<dimensionality>::Coordinate_>, std::placeholders::_1, std::placeholders::_2, tolerance));
}

template<int dimensionality>
bool operator==(VectorSpace<dimensionality> const &lhs, VectorSpace<dimensionality> const &rhs) {
  return IsEqual(lhs, rhs);
}

template<int dimensionality>
typename VectorSpace<dimensionality>::Coordinate_ const &
VectorSpace<dimensionality>::operator[](Index const &coordinate) const {
#ifndef NDEBUG
  try {
    ThrowIfIndexIsInvalid(coordinate);
  } catch (OutOfRange const &exception) {
    Throw(exception, "splinelib::sources::vector_spaces::VectorSpace::operator[]");
  }
#endif
  return coordinates_[coordinate.Get()];
}

template<int dimensionality>
int VectorSpace<dimensionality>::GetNumberOfCoordinates() const {
  return coordinates_.size();
}

template<int dimensionality>
void VectorSpace<dimensionality>::Replace(Index const &coordinate_index, Coordinate_ coordinate) {
#ifndef NDEBUG
  try {
    ThrowIfIndexIsInvalid(coordinate_index);
  } catch (OutOfRange const &exception) { Throw(exception, "splinelib::sources::vector_spaces::VectorSpace::Replace"); }
#endif
  coordinates_[coordinate_index.Get()] = std::move(coordinate);
}

template<int dimensionality>
void VectorSpace<dimensionality>::Insert(Index const &coordinate_index, Coordinate_ coordinate) {
#ifndef NDEBUG
  try {
    ThrowIfIndexIsInvalid(coordinate_index);
  } catch (OutOfRange const &exception) { Throw(exception, "splinelib::sources::vector_spaces::VectorSpace::Insert"); }
#endif
  coordinates_.insert(coordinates_.begin() + coordinate_index.Get(), std::move(coordinate));
}

template<int dimensionality>
void VectorSpace<dimensionality>::Erase(Index const &coordinate_index) {
#ifndef NDEBUG
  try {
    ThrowIfIndexIsInvalid(coordinate_index);
  } catch (OutOfRange const &exception) { Throw(exception, "splinelib::sources::vector_spaces::VectorSpace::Erase"); }
#endif
  coordinates_.erase(coordinates_.begin() + coordinate_index.Get());
}

template<int dimensionality>
Coordinate VectorSpace<dimensionality>::DetermineMaximumDistanceFromOrigin(Tolerance const &tolerance) const {
#ifndef NDEBUG
  try {
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const &exception) {
    Throw(exception, "splinelib::sources::vector_spaces::VectorSpace::DetermineMaximumDistanceFromOrigin");
  }
#endif
  Coordinate maximum_distance{};
  std::for_each(coordinates_.begin(), coordinates_.end(), [&] (Coordinate_ const &coordinate) {
      maximum_distance = std::max(utilities::std_container_operations::TwoNorm(coordinate), maximum_distance,
                                  std::bind(utilities::IsLess<typename Coordinate::Name_, typename Coordinate::Type_>,
                                            std::placeholders::_1, std::placeholders::_2, tolerance)); });
  return maximum_distance;
}

template<int dimensionality>
typename VectorSpace<dimensionality>::OutputInformation_ VectorSpace<dimensionality>::Write(Precision const &precision)
    const {
  return OutputInformation_{utilities::string_operations::Write<std::tuple_element_t<0, OutputInformation_>>(
                                                              coordinates_, precision)};
}

#ifndef NDEBUG
template<int dimensionality>
void VectorSpace<dimensionality>::ThrowIfIndexIsInvalid(Index const &coordinate) const {
  Index::ThrowIfNamedIntegerIsOutOfBounds(coordinate, coordinates_.size() - 1);
}
#endif
