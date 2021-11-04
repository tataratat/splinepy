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

#ifndef SOURCES_SPLINES_SPLINE_ITEM_HPP_
#define SOURCES_SPLINES_SPLINE_ITEM_HPP_

#include "Sources/Utilities/named_type.hpp"

namespace splinelib::sources::splines {

// SplineItems allow to store lists of splines without losing the information on their parametric dimensionality,
// dimensionality or whether they are rational or not.
//
// Example:
//   using std::make_shared;
//   Vector<SplineItem> const splines{make_shared<BSpline<2, 3>>(), make_shared<Nurbs<3, 4>>()};
//   int const &two = splines[0]->parametric_dimensionality_;
class SplineItem {
 public:
  using Type_ = Type;

  virtual ~SplineItem() = default;

  friend bool IsEqual(SplineItem const &lhs, SplineItem const &rhs, Tolerance const &tolerance);
  friend bool operator==(SplineItem const &lhs, SplineItem const &rhs);

  int dimensionality_, parametric_dimensionality_;
  bool is_rational_;

 protected:
  SplineItem() = default;
  SplineItem(int parametric_dimensionality, int dimensionality, bool is_rational);
  SplineItem(SplineItem const &other) = default;
  SplineItem(SplineItem &&other) noexcept = default;
  SplineItem & operator=(SplineItem const &rhs) = default;
  SplineItem & operator=(SplineItem &&rhs) noexcept = default;
};

bool IsEqual(SplineItem const &lhs, SplineItem const &rhs, Tolerance const &tolerance = kEpsilon);
bool operator==(SplineItem const &lhs, SplineItem const &rhs);

}  // namespace splinelib::sources::splines

#endif  // SOURCES_SPLINES_SPLINE_ITEM_HPP_
