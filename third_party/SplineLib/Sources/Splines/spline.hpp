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

#ifndef SOURCES_SPLINES_SPLINE_HPP_
#define SOURCES_SPLINES_SPLINE_HPP_

#include <algorithm>
#include <utility>

#include "Sources/ParameterSpaces/parameter_space.hpp"
#include "Sources/Splines/spline_item.hpp"
#include "Sources/Utilities/error_handling.hpp"
#include "Sources/Utilities/named_type.hpp"
#include "Sources/Utilities/numeric_operations.hpp"
#include "Sources/Utilities/std_container_operations.hpp"
#include "Sources/VectorSpaces/vector_space.hpp"

namespace splinelib::sources::splines {

template<int parametric_dimensionality, int dimensionality> class Spline;

template<int parametric_dimensionality, int dimensionality>
bool IsEqual(Spline<parametric_dimensionality, dimensionality> const &lhs,
             Spline<parametric_dimensionality, dimensionality> const &rhs, Tolerance const &tolerance = kEpsilon);
template<int parametric_dimensionality, int dimensionality>
bool operator==(Spline<parametric_dimensionality, dimensionality> const &lhs,
                Spline<parametric_dimensionality, dimensionality> const &rhs);

// Splines are (non-)rational mappings from parameter spaces of arbitrary parametric_dimensionality to vector spaces of
// arbitrary dimensionality.  They can be used, e.g., to seamlessly integrate geometry representation, computational
// analysis, and optimization.
//
// Example (see, e.g., NURBS book Exe. 3.8 or Exe. 4.4):
//   spline.RefineKnots(Dimension{1}, {Spline<2, 3>::Knot_{0.5}});
template<int parametric_dimensionality, int dimensionality>
class Spline : public SplineItem {
 protected:
  using VectorSpace_ = vector_spaces::VectorSpace<dimensionality>;

 public:
  using Base_ = SplineItem;
  using Coordinate_ = typename VectorSpace_::Coordinate_;
  using Coordinates_ = typename VectorSpace_::Coordinates_;
  using ParameterSpace_ = parameter_spaces::ParameterSpace<parametric_dimensionality>;
  using Derivative_ = typename ParameterSpace_::Derivative_;
  using Knot_ = typename ParameterSpace_::Knot_;
  using Knots_ = typename ParameterSpace_::Knots_;
  using NumberOfParametricCoordinates_ = typename ParameterSpace_::NumberOfParametricCoordinates_;
  using ParametricCoordinate_ = typename ParameterSpace_::ParametricCoordinate_;

  ~Spline() override = default;

  // Comparison based on tolerance.
  friend bool IsEqual<parametric_dimensionality, dimensionality>(Spline const &lhs, Spline const &rhs,
                                                                 Tolerance const &tolerance);
  // Comparison based on numeric_operations::GetEpsilon<Tolerance>().
  friend bool operator==<parametric_dimensionality, dimensionality>(Spline const &lhs, Spline const &rhs);
  virtual Coordinate_ operator()(ParametricCoordinate_ const &parametric_coordinate,
                                 Tolerance const &tolerance = kEpsilon) const = 0;
  virtual Coordinate_ operator()(ParametricCoordinate_ const &parametric_coordinate, Derivative_ const &derivative,
                                 Tolerance const &tolerance = kEpsilon) const = 0;

  virtual void InsertKnot(Dimension const &dimension, Knot_ knot, Multiplicity const &multiplicity = kMultiplicity,
                          Tolerance const &tolerance = kEpsilon) const = 0;
  void RefineKnots(Dimension const &dimension, Knots_ knots, Multiplicity const &multiplicity = kMultiplicity,
                   Tolerance const &tolerance = kEpsilon) const;
  virtual Multiplicity RemoveKnot(Dimension const &dimension, Knot_ const &knot, Tolerance const &tolerance_removal,
      Multiplicity const &multiplicity = kMultiplicity, Tolerance const &tolerance = kEpsilon) const = 0;
  Multiplicity CoarsenKnots(Dimension const &dimension, Knots_ const &knots, Tolerance const &tolerance_removal,
      Multiplicity const &multiplicity = kMultiplicity, Tolerance const &tolerance = kEpsilon) const;
  virtual void ElevateDegree(Dimension const &dimension, Multiplicity const &multiplicity = kMultiplicity,
                             Tolerance const &tolerance = kEpsilon) const = 0;
  virtual bool ReduceDegree(Dimension const &dimension, Tolerance const &tolerance_reduction,
      Multiplicity const &multiplicity = kMultiplicity, Tolerance const &tolerance = kEpsilon) const = 0;

  virtual Coordinate ComputeUpperBoundForMaximumDistanceFromOrigin(Tolerance const &tolerance = kEpsilon) const = 0;
  virtual Coordinates_ Sample(NumberOfParametricCoordinates_ const &number_of_parametric_coordinates,
                              Tolerance const &tolerance = kEpsilon) const;

 protected:
  using Index_ = typename ParameterSpace_::Index_;

  Spline() = default;
  explicit Spline(bool is_rational);
  Spline(SharedPointer<ParameterSpace_> parameter_space, bool is_rational);
  Spline(Spline const &other);
  Spline(Spline &&other) noexcept = default;
  Spline & operator=(Spline const &rhs);
  Spline & operator=(Spline &&rhs) noexcept = default;

  SharedPointer<ParameterSpace_> parameter_space_;
};

#include "Sources/Splines/spline.inc"

}  // namespace splinelib::sources::splines

#endif  // SOURCES_SPLINES_SPLINE_HPP_
