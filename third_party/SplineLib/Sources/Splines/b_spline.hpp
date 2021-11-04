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

#ifndef SOURCES_SPLINES_B_SPLINE_HPP_
#define SOURCES_SPLINES_B_SPLINE_HPP_

#include <algorithm>
#include <iterator>
#include <utility>

#include "Sources/Splines/spline.hpp"
#include "Sources/Utilities/error_handling.hpp"
#include "Sources/Utilities/index.hpp"
#include "Sources/Utilities/named_type.hpp"
#include "Sources/Utilities/numeric_operations.hpp"
#include "Sources/Utilities/std_container_operations.hpp"
#include "Sources/VectorSpaces/vector_space.hpp"

namespace splinelib::sources::splines {

template<int parametric_dimensionality, int dimensionality> class BSpline;

template<int parametric_dimensionality, int dimensionality>
bool IsEqual(BSpline<parametric_dimensionality, dimensionality> const &lhs,
             BSpline<parametric_dimensionality, dimensionality> const &rhs, Tolerance const &tolerance = kEpsilon);
template<int parametric_dimensionality, int dimensionality>
bool operator==(BSpline<parametric_dimensionality, dimensionality> const &lhs,
                BSpline<parametric_dimensionality, dimensionality> const &rhs);

// B-splines are non-rational splines.  Currently only single-patch B-splines are supported.
//
// Example (see, e.g., NURBS book Exe. 3.8):
//   using Surface = BSpline<2, 3>;
//   Surface::Coordinate_ const &coordinate0 = surface({});  // Evaluating the B-spline S(0.0, 0.0) results in P_{0,0}.
//   Multiplicity const &successful_removals = surface.RemoveKnot(Dimension{}, Surface::Knot_{0.5});
//   bool const &successful = surface.ReduceDegree(Dimension{1}, kEpsilon);  // True if spline's degree p_0 be reduced.
template<int parametric_dimensionality, int dimensionality>
class BSpline : public Spline<parametric_dimensionality, dimensionality> {
 public:
  using Base_ = Spline<parametric_dimensionality, dimensionality>;
  using Coordinate_ = typename Base_::Coordinate_;
  using Derivative_ = typename Base_::Derivative_;
  using Knot_ = typename Base_::Knot_;
  using ParameterSpace_ = typename Base_::ParameterSpace_;
  using ParametricCoordinate_ = typename Base_::ParametricCoordinate_;
  using VectorSpace_ = typename Base_::VectorSpace_;
  using OutputInformation_ = Tuple<typename ParameterSpace_::OutputInformation_,
                                   typename VectorSpace_::OutputInformation_>;

  BSpline();
  BSpline(SharedPointer<ParameterSpace_> parameter_space, SharedPointer<VectorSpace_> vector_space);
  BSpline(BSpline const &other);
  BSpline(BSpline &&other) noexcept = default;
  BSpline & operator=(BSpline const &rhs);
  BSpline & operator=(BSpline &&rhs) noexcept = default;
  ~BSpline() override = default;

  // Comparison based on tolerance.
  friend bool IsEqual<parametric_dimensionality, dimensionality>(BSpline const &lhs, BSpline const &rhs,
                                                                 Tolerance const &tolerance);
  // Comparison based on numeric_operations::GetEpsilon<Tolerance>().
  friend bool operator==<parametric_dimensionality, dimensionality>(BSpline const &lhs, BSpline const &rhs);
  Coordinate_ operator()(ParametricCoordinate_ const &parametric_coordinate, Tolerance const &tolerance = kEpsilon)
      const override;
  Coordinate_ operator()(ParametricCoordinate_ const &parametric_coordinate, Derivative_ const &derivative,
                         Tolerance const &tolerance = kEpsilon) const override;

  void InsertKnot(Dimension const &dimension, Knot_ knot, Multiplicity const &multiplicity = kMultiplicity,
                  Tolerance const &tolerance = kEpsilon) const override;
  // Tries to interpret knot removal as the inverse process of knot insertion.
  Multiplicity RemoveKnot(Dimension const &dimension, Knot_ const &knot, Tolerance const &tolerance_removal,
      Multiplicity const &multiplicity = kMultiplicity, Tolerance const &tolerance = kEpsilon) const override;
  void ElevateDegree(Dimension const &dimension, Multiplicity const &multiplicity = kMultiplicity,
                     Tolerance const &tolerance = kEpsilon) const override;
  // Tries to interpret degree reduction as the inverse process of degree elevation.
  bool ReduceDegree(Dimension const &dimension, Tolerance const &tolerance_reduction,
      Multiplicity const &multiplicity = kMultiplicity, Tolerance const &tolerance = kEpsilon) const override;

  Coordinate ComputeUpperBoundForMaximumDistanceFromOrigin(Tolerance const &tolerance = kEpsilon) const override;
  OutputInformation_ Write(Precision const &precision = kPrecision) const;

 protected:
  SharedPointer<VectorSpace_> vector_space_;

 private:
  using BezierInformation_ = typename ParameterSpace_::BezierInformation_;
  using BinomialRatios_ = typename ParameterSpace_::BinomialRatios_;
  using Index_ = typename Base_::Index_;
  using IndexLength_ = typename Index_::Length_;
  using IndexValue_ = typename Index_::Value_;
  using KnotRatios_ = typename ParameterSpace_::KnotRatios_;
  using Knots_ = typename Base_::Knots_;
  using BinomialRatio_ = typename BinomialRatios_::value_type;
  using KnotRatio_ = typename KnotRatios_::value_type;

  BezierInformation_ MakeBezier(Dimension const &dimension, Tolerance const &tolerance = kEpsilon) const;
};

#include "Sources/Splines/b_spline.inc"

}  // namespace splinelib::sources::splines

#endif  // SOURCES_SPLINES_B_SPLINE_HPP_
