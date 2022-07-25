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

#ifndef SOURCES_PARAMETERSPACES_KNOT_VECTOR_HPP_
#define SOURCES_PARAMETERSPACES_KNOT_VECTOR_HPP_

#include "Sources/Utilities/named_type.hpp"
#include "Sources/Utilities/std_container_operations.hpp"
#include "Sources/Utilities/string_operations.hpp"

namespace splinelib::sources::parameter_spaces {

// KnotVectors are sequences of non-decreasing real numbers (called knots).
//
// Example:
//   using Knot = KnotVector::Knot_;
//   constexpr Knot const kStart{}, kEnd{1.0};
//   KnotVector knot_vector{kStart, kStart, Knot{0.5}, kEnd, kEnd};  // Knot vector U is {0.0, 0.0, 0.5, 1.0, 1.0}.
//   int const &number_of_knots = knot_vector.GetSize();  // Stores m+1 = 5 in number_of_knots.
//   Knot const &zero_point_five = knot_vector[2];  // Third knot u is 0.5.
//   ParametricCoordinate const parametric_coordinate{0.25};
//   KnotSpan const &knot_span_0_25 = knot_vector.FindSpan(parametric_coordinate),  // Knot span i = 1.
//                  &knot_span_1_0 = knot_vector.FindSpan(knot);  // The last knot span i = 3 is treated special.
//   Multiplicity const &zero = knot_vector.DetermineMultiplicity(parametric_coordinate);  // Multiplicity s = 0.
//   knot_vector.DoesParametricCoordinateEqualBack(knot);  // Evaluates to true as 1.0 is the last knot.
class KnotVector {
 public:
  using Knot_ = ParametricCoordinate;
  using OutputInformation_ = StringVector;
  using Knots_ = Vector<Knot_>;
  using Type_ = Knot_::Type_;

  KnotVector() = default;
  explicit KnotVector(Knots_ knots, Tolerance const &tolerance = kEpsilon);
  KnotVector(KnotVector const &other) = default;
  KnotVector(KnotVector &&other) noexcept = default;
  KnotVector & operator=(KnotVector const &rhs) = default;
  KnotVector & operator=(KnotVector &&rhs) noexcept = default;
  virtual ~KnotVector() = default;

  // Comparison based on tolerance.
  friend bool IsEqual(KnotVector const &lhs, KnotVector const &rhs, Tolerance const &tolerance);
  // Comparison based on numeric_operations::GetEpsilon<Tolerance>().
  friend bool operator==(KnotVector const &lhs, KnotVector const &rhs);
  virtual Knot_ const & operator[](Index const &index) const;

  virtual int GetSize() const;
  virtual Knot_ const & GetFront() const;
  virtual Knot_ const & GetBack() const;

  virtual bool DoesParametricCoordinateEqualBack(ParametricCoordinate const &parametric_coordinate,
                                                 Tolerance const &tolerance = kEpsilon) const;
  virtual bool DoesParametricCoordinateEqualFrontOrBack(ParametricCoordinate const &parametric_coordinate,
                                                        Tolerance const &tolerance = kEpsilon) const;
  // Returns the (m-s)th or ith knot span if the parametric coordinate u equals the last knot u_m of multiplicity s or
  // is in the interval [u_i, u_{i+1}) (implying the knot span is non-zero, i.e., u_i < u_{i+1}), respectively.
  virtual KnotSpan FindSpan(ParametricCoordinate const &parametric_coordinate, Tolerance const &tolerance = kEpsilon)
      const;
  virtual Multiplicity DetermineMultiplicity(Knot_ const &knot, Tolerance const &tolerance = kEpsilon) const;
  virtual Knots_ GetUniqueKnots(Tolerance const &tolerance = kEpsilon) const;

  virtual void Insert(Knot_ knot, Multiplicity const &multiplicity = kMultiplicity,
                      Tolerance const &tolerance = kEpsilon);
  virtual void Remove(Knot_ const &knot, Multiplicity const &multiplicity = kMultiplicity,
                      Tolerance const &tolerance = kEpsilon);
  virtual void IncreaseMultiplicities(Multiplicity const &multiplicity = kMultiplicity,
                                      Tolerance const &tolerance = kEpsilon);
  virtual void DecreaseMultiplicities(Multiplicity const &multiplicity = kMultiplicity,
                                      Tolerance const &tolerance = kEpsilon);

  virtual OutputInformation_ Write(Precision const &precision = kPrecision) const;

#ifndef NDEBUG
  virtual void ThrowIfParametricCoordinateIsOutsideScope(ParametricCoordinate const &parametric_coordinate,
                                                         Tolerance const &tolerance = kEpsilon) const;
#endif

 protected:
  Knots_ knots_;

 private:
  using ConstIterator_ = typename Knots_::const_iterator;

#ifndef NDEBUG
  void ThrowIfTooSmallOrNotNonDecreasing(Tolerance const &tolerance) const;
#endif
};

bool IsEqual(KnotVector const &lhs, KnotVector const &rhs, Tolerance const &tolerance = kEpsilon);
bool operator==(KnotVector const &lhs, KnotVector const &rhs);

template<int parametric_dimensionality>
using KnotVectors = Array<SharedPointer<KnotVector>, parametric_dimensionality>;
template<int parametric_dimensionality>
using KnotVectorsOutputInformation = Array<KnotVector::OutputInformation_, parametric_dimensionality>;
}  // namespace splinelib::sources::parameter_spaces

#endif  // SOURCES_PARAMETERSPACES_KNOT_VECTOR_HPP_
