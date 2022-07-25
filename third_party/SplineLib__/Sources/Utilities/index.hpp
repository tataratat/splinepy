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

#ifndef SOURCES_UTILITIES_INDEX_HPP_
#define SOURCES_UTILITIES_INDEX_HPP_

#include <algorithm>
#include <functional>
#include <numeric>
#include <utility>

#include "Sources/Utilities/error_handling.hpp"
#include "Sources/Utilities/named_type.hpp"
#include "Sources/Utilities/std_container_operations.hpp"

namespace splinelib::sources::utilities {

template<int size> class Index;

template<int size>
Index<size> operator+(Index<size> const &lhs, typename Index<size>::Value_ const &rhs);
template<int size>
Index<size> operator-(Index<size> const &lhs, typename Index<size>::Value_ const &rhs);
template<int size>
bool operator==(Index<size> const &lhs, Index<size> const &rhs);
template<int size>
bool operator!=(Index<size> const &lhs, Index<size> const &rhs);

// Indices map (column-major order) size-dimensional arrays to one-dimensional arrays by mapping indices of size to
// single integers (and vice versa).
//
// Example:
//   using TripleIndex = Index<3>;
//   TripleIndex::Length_ length{Length{3}, Length{2}, Length{4}};
//   TripleIndex index{length};  // New valid index with initial value of {0, 0, 0}.
//   int const &twenty_four = index.GetTotalNumberOfIndices();
//   ++index;  // Moves index from value {0, 0, 0} to value {1, 0, 0}.
//   using Index1d = splinelib::Index;
//   constexpr Index1d const k1{1};
//   index += TripleIndex::Value_{k1, k1, Index1d{3}};  // Moves index to its maximum {2, 1, 3}.
//   Index1d const &twenty_three = index.GetIndex1d();  // The maximum one-dimensional index is 23.
//   ++index;  // Index becomes invalid when moving from {2, 1, 3} to {0, 0, 0}.
//   index == TripleIndex::Behind(length});  // Evaluates to true as the index is behind.
template<int size>
class Index {
 private:
  template<typename Type>
  using Array_ = Array<Type, size>;

 public:
  using Index_ = splinelib::Index;
  using Length_ = Array_<Length>;
  using Value_ = Array_<Index_>;

  Index() = default;
  explicit Index(Length_ length, Value_ value = Value_{}, bool invalid = false);
  Index(Index const &other) = default;
  Index(Index &&other) noexcept = default;
  Index & operator=(Index const &rhs) = default;
  Index & operator=(Index &&rhs) noexcept = default;
  virtual ~Index() = default;

  Index & operator+=(Value_ const &value);
  Index & operator-=(Value_ const &value);
  Index & operator++();
  Index & Increment(Dimension const &dimension);
  Index & operator--();
  Index & Decrement(Dimension const &dimension);
  friend Index operator+<size>(Index const &lhs, Value_ const &rhs);
  friend Index operator-<size>(Index const &lhs, Value_ const &rhs);
  friend bool operator==<size>(Index const &lhs, Index const &rhs);
  friend bool operator!=<size>(Index const &lhs, Index const &rhs);
  virtual Index_ const & operator[](Dimension const &dimension) const;
  static Index First(Length_ length);
  static Index Behind(Length_ length);
  static Index Last(Length_ length);
  static Index Before(Length_ length);

  int GetTotalNumberOfIndices() const;
  Value_ GetIndex() const;
  Index_ GetIndex1d() const;

 protected:
  bool invalid_;
  Length_ length_;
  Value_ value_;

 private:
  int DetermineStride(Length_ const &length, Dimension const &dimension) const;

#ifndef NDEBUG
  static void ThrowIfValueIsInvalid(Length_ const &length, Value_ const &value);
#endif
};

#include "Sources/Utilities/index.inc"

}  // namespace splinelib::sources::utilities

#endif  // SOURCES_UTILITIES_INDEX_HPP_
