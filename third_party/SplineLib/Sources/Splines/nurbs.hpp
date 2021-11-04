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

#ifndef SOURCES_SPLINES_NURBS_HPP_
#define SOURCES_SPLINES_NURBS_HPP_

#include <algorithm>
#include <deque>
#include <iterator>
#include <utility>

#include "Sources/ParameterSpaces/parameter_space.hpp"
#include "Sources/Splines/b_spline.hpp"
#include "Sources/Splines/spline.hpp"
#include "Sources/Utilities/error_handling.hpp"
#include "Sources/Utilities/math_operations.hpp"
#include "Sources/Utilities/named_type.hpp"
#include "Sources/Utilities/numeric_operations.hpp"
#include "Sources/Utilities/std_container_operations.hpp"
#include "Sources/VectorSpaces/weighted_vector_space.hpp"

namespace splinelib::sources::splines {

template<int parametric_dimensionality, int dimensionality> class Nurbs;

template<int parametric_dimensionality, int dimensionality>
bool IsEqual(Nurbs<parametric_dimensionality, dimensionality> const &lhs,
             Nurbs<parametric_dimensionality, dimensionality> const &rhs, Tolerance const &tolerance = kEpsilon);
template<int parametric_dimensionality, int dimensionality>
bool operator==(Nurbs<parametric_dimensionality, dimensionality> const &lhs,
                Nurbs<parametric_dimensionality, dimensionality> const &rhs);

// NURBSs are rational B-splines.  Currently only single-patch NURBSs are supported.
//
// Example (see NURBS book Exe. 4.4):
//   using Curve = Nurbs<1, 3>;
//   Curve::Coordinate_ const &first_derivative0 = curve({}, Derivative{1});  // Evaluating S'(0.0, 0.0) (Eq. (4.9)).
//   curve.InsertKnot(Dimension{1}, Curve::Knot_{0.5});  // Insert u = 0.5 into U_1.
//   curve.ElevateDegree(Dimension{});  // Raise the spline's degree p_0 by one.
template<int parametric_dimensionality, int dimensionality>
class Nurbs : public Spline<parametric_dimensionality, dimensionality> {
 public:
  using Base_ = Spline<parametric_dimensionality, dimensionality>;
  using Coordinate_ = typename Base_::Coordinate_;
  using Derivative_ = typename Base_::Derivative_;
  using Knot_ = typename Base_::Knot_;
  using ParameterSpace_ = typename Base_::ParameterSpace_;
  using ParametricCoordinate_ = typename Base_::ParametricCoordinate_;
  using WeightedVectorSpace_ = vector_spaces::WeightedVectorSpace<dimensionality>;
  using OutputInformation_ = Tuple<typename ParameterSpace_::OutputInformation_,
                                   typename WeightedVectorSpace_::OutputInformation_>;

  Nurbs();
  Nurbs(SharedPointer<ParameterSpace_> parameter_space, SharedPointer<WeightedVectorSpace_> weighted_vector_space);
  Nurbs(Nurbs const &other);
  Nurbs(Nurbs &&other) noexcept = default;
  Nurbs & operator=(Nurbs const &rhs);
  Nurbs & operator=(Nurbs &&rhs) noexcept = default;
  ~Nurbs() override = default;

  // Comparison based on tolerance.
  friend bool IsEqual<parametric_dimensionality, dimensionality>(Nurbs const &lhs, Nurbs const &rhs,
                                                                 Tolerance const &tolerance);
  // Comparison based on numeric_operations::GetEpsilon<Tolerance>().
  friend bool operator==<parametric_dimensionality, dimensionality>(Nurbs const &lhs, Nurbs const &rhs);
  Coordinate_ operator()(ParametricCoordinate_ const &parametric_coordinate,
                         Tolerance const &tolerance = kEpsilon) const final;
  Coordinate_ operator()(ParametricCoordinate_ const &parametric_coordinate, Derivative_ const &derivative,
                         Tolerance const &tolerance = kEpsilon) const final;

  void InsertKnot(Dimension const &dimension, Knot_ knot, Multiplicity const &multiplicity = kMultiplicity,
                  Tolerance const &tolerance = kEpsilon) const final;
  Multiplicity RemoveKnot(Dimension const &dimension, Knot_ const &knot, Tolerance const &tolerance_removal,
      Multiplicity const &multiplicity = kMultiplicity, Tolerance const &tolerance = kEpsilon) const final;
  void ElevateDegree(Dimension const &dimension, Multiplicity const &multiplicity = kMultiplicity,
                     Tolerance const &tolerance = kEpsilon) const final;
  bool ReduceDegree(Dimension const &dimension, Tolerance const &tolerance_removal,
      Multiplicity const &multiplicity = kMultiplicity, Tolerance const &tolerance = kEpsilon) const final;

  Coordinate ComputeUpperBoundForMaximumDistanceFromOrigin(Tolerance const &tolerance = kEpsilon) const final;
  OutputInformation_ Write(Precision const &precision = kPrecision) const;

 protected:
  using HomogeneousBSpline_ = BSpline<parametric_dimensionality, dimensionality + 1>;

  SharedPointer<HomogeneousBSpline_> homogeneous_b_spline_;
  SharedPointer<WeightedVectorSpace_> weighted_vector_space_;
};

#include "Sources/Splines/nurbs.inc"

}  // namespace splinelib::sources::splines

#endif  // SOURCES_SPLINES_NURBS_HPP_
