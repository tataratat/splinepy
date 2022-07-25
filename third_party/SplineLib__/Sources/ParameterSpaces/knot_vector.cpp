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

#include "Sources/ParameterSpaces/knot_vector.hpp"

#include <algorithm>
#include <functional>
#include <iterator>
#include <utility>

#include "Sources/Utilities/error_handling.hpp"
#include "Sources/Utilities/numeric_operations.hpp"
#include "Sources/Utilities/std_container_operations.hpp"

namespace splinelib::sources::parameter_spaces {

using Knot = KnotVector::Knot_;
using std::move, std::to_string;
#ifndef NDEBUG
using utilities::numeric_operations::ThrowIfToleranceIsNegative;
#endif

KnotVector::KnotVector(Knots_ knots, Tolerance const &tolerance) : knots_(move(knots)) {
#ifndef NDEBUG
  Message const kName{"splinelib::sources::parameter_spaces::KnotVector::KnotVector"};

  try {
    ThrowIfToleranceIsNegative(tolerance);
    ThrowIfTooSmallOrNotNonDecreasing(tolerance);
  } catch (DomainError const &exception) { Throw(exception, kName); }
    catch (InvalidArgument const &exception) { Throw(exception, kName); }
#endif
}

bool IsEqual(KnotVector const &lhs, KnotVector const &rhs, Tolerance const &tolerance) {
#ifndef NDEBUG
  try {
    ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const &exception) {
    Throw(exception, "splinelib::sources::parameter_spaces::IsEqual::KnotVector");
  }
#endif
  return utilities::std_container_operations::DoesContainEqualValues(lhs.knots_, rhs.knots_, tolerance);
}

bool operator==(KnotVector const &lhs, KnotVector const &rhs) {
  return IsEqual(lhs, rhs);
}

Knot const & KnotVector::operator[](Index const &index) const {
#ifndef NDEBUG
  try {
    Index::ThrowIfNamedIntegerIsOutOfBounds(index, knots_.size() - 1);
  } catch (OutOfRange const &exception) {
    Throw(exception, "splinelib::sources::parameter_spaces::KnotVector::operator[]");
  }
#endif
  return knots_[index.Get()];
}

int KnotVector::GetSize() const {
  return knots_.size();
}

Knot const & KnotVector::GetFront() const {
  return knots_[0];
}

Knot const & KnotVector::GetBack() const {
  return knots_[knots_.size() - 1];
}

bool KnotVector::DoesParametricCoordinateEqualBack(ParametricCoordinate const &parametric_coordinate,
                                                   Tolerance const &tolerance) const {
#ifndef NDEBUG
  Message const kName{"splinelib::sources::parameter_spaces::KnotVector::DoesParametricCoordinateEqualBack"};

  try {
    ThrowIfToleranceIsNegative(tolerance);
    ThrowIfTooSmallOrNotNonDecreasing(tolerance);
  } catch (DomainError const &exception) { Throw(exception, kName); }
    catch (InvalidArgument const &exception) { Throw(exception, kName); }
#endif
  return IsEqual(parametric_coordinate, GetBack(), tolerance);
}

bool KnotVector::DoesParametricCoordinateEqualFrontOrBack(ParametricCoordinate const &parametric_coordinate,
                                                          Tolerance const &tolerance) const {
#ifndef NDEBUG
  Message const kName{"splinelib::sources::parameter_spaces::KnotVector::DoesParametricCoordinateEqualFrontOrBack"};

  try {
    ThrowIfToleranceIsNegative(tolerance);
    ThrowIfTooSmallOrNotNonDecreasing(tolerance);
  } catch (DomainError const &exception) { Throw(exception, kName); }
    catch (InvalidArgument const &exception) { Throw(exception, kName); }
#endif
  return (IsEqual(parametric_coordinate, GetFront(), tolerance) ||
          DoesParametricCoordinateEqualBack(parametric_coordinate, tolerance));
}

KnotSpan KnotVector::FindSpan(ParametricCoordinate const &parametric_coordinate, Tolerance const &tolerance) const {
#ifndef NDEBUG
  Message const kName{"splinelib::sources::parameter_spaces::KnotVector::FindSpan"};

  try {
    ThrowIfToleranceIsNegative(tolerance);
    ThrowIfParametricCoordinateIsOutsideScope(parametric_coordinate, tolerance);
    ThrowIfTooSmallOrNotNonDecreasing(tolerance);
  } catch (DomainError const &exception) { Throw(exception, kName); }
    catch (InvalidArgument const &exception) { Throw(exception, kName); }
    catch (OutOfRange const &exception) { Throw(exception, kName); }
#endif
  ConstIterator_ const &knots_begin = knots_.begin(), &knots_end = knots_.end();
  std::function<bool(ParametricCoordinate, Knot_)> const &is_less =
      [&] (Knot_ const &knot, ParametricCoordinate const &parametric_coordinate) {
          return IsLess(knot, parametric_coordinate, tolerance); };
  return KnotSpan{static_cast<int>(std::distance(knots_begin, DoesParametricCoordinateEqualBack(parametric_coordinate,
                      tolerance) ? std::lower_bound(knots_begin, knots_end, parametric_coordinate, is_less) :
                                   std::upper_bound(knots_begin, knots_end, parametric_coordinate, is_less)) - 1)};
}

Multiplicity KnotVector::DetermineMultiplicity(ParametricCoordinate const &parametric_coordinate,
                                               Tolerance const &tolerance) const {
#ifndef NDEBUG
  Message const kName{"splinelib::sources::parameter_spaces::KnotVector::DetermineMultiplicity"};

  try {
    ThrowIfToleranceIsNegative(tolerance);
    ThrowIfTooSmallOrNotNonDecreasing(tolerance);
  } catch (DomainError const &exception) { Throw(exception, kName); }
    catch (InvalidArgument const &exception) { Throw(exception, kName); }
#endif
  return Multiplicity{static_cast<int>(std::count_if(knots_.begin(), knots_.end(), [&] (Knot_ const &current_knot) {
                                           return IsEqual(current_knot, parametric_coordinate, tolerance); }))};
}

KnotVector::Knots_ KnotVector::GetUniqueKnots(Tolerance const &tolerance) const {
#ifndef NDEBUG
  try {
    ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const &exception) {
    Throw(exception, "splinelib::sources::parameter_spaces::KnotVector::GetUniqueKnots");
  }
#endif
  Knots_ unique_knots;
  std::unique_copy(knots_.begin(), knots_.end(), std::back_inserter(unique_knots),
      [&] (Knot_ const &lhs_knot, Knot_ const &rhs_knot) { return IsEqual(lhs_knot, rhs_knot, tolerance); });
  return unique_knots;
}

void KnotVector::Insert(Knot_ knot, Multiplicity const &multiplicity, Tolerance const &tolerance) {
#ifndef NDEBUG
  Message const kName{"splinelib::sources::parameter_spaces::KnotVector::Insert"};

  try {
    ThrowIfToleranceIsNegative(tolerance);
    ThrowIfParametricCoordinateIsOutsideScope(knot, tolerance);
    ThrowIfTooSmallOrNotNonDecreasing(tolerance);
  } catch (DomainError const &exception) { Throw(exception, kName); }
    catch (InvalidArgument const &exception) { Throw(exception, kName); }
    catch (OutOfRange const &exception) { Throw(exception, kName); }
#endif
  knots_.insert(knots_.begin() + FindSpan(knot, tolerance).Get() + 1, multiplicity.Get(), move(knot));
}

void KnotVector::Remove(Knot_ const &knot, Multiplicity const &multiplicity, Tolerance const &tolerance) {
#ifndef NDEBUG
  Message const kName{"splinelib::sources::parameter_spaces::KnotVector::Remove"};

  try {
    ThrowIfToleranceIsNegative(tolerance);
    ThrowIfTooSmallOrNotNonDecreasing(tolerance);
  } catch (DomainError const &exception) { Throw(exception, kName); }
    catch (InvalidArgument const &exception) { Throw(exception, kName); }
#endif
  Multiplicity::Type_ const number_of_removals{std::min(multiplicity.Get(), DetermineMultiplicity(knot,
                                                                                                  tolerance).Get())};
  if (number_of_removals != 0) {
    KnotSpan const &knot_span = FindSpan(knot, tolerance);
    if (DoesParametricCoordinateEqualBack(knot, tolerance)) {
      ConstIterator_ const &end = knots_.end();
      knots_.erase(end - number_of_removals, end);
    } else {
      ConstIterator_ const &first_knot = (knots_.begin() + knot_span.Get());
      knots_.erase(first_knot - (number_of_removals - 1), first_knot + 1);
    }
  }
}

void KnotVector::IncreaseMultiplicities(Multiplicity const &multiplicity, Tolerance const &tolerance) {
  for (Knot_ const &knot : GetUniqueKnots(tolerance)) Insert(knot, multiplicity, tolerance);
}

void KnotVector::DecreaseMultiplicities(Multiplicity const &multiplicity, Tolerance const &tolerance) {
  if (GetSize() > 2) for (Knot_ const &knot : GetUniqueKnots(tolerance)) Remove(knot, multiplicity, tolerance);
}

typename KnotVector::OutputInformation_ KnotVector::Write(Precision const &precision) const {
  return utilities::string_operations::Write<OutputInformation_>(knots_, precision);
}

#ifndef NDEBUG
void KnotVector::ThrowIfParametricCoordinateIsOutsideScope(ParametricCoordinate const &parametric_coordinate,
                                                           Tolerance const &tolerance) const {
  Knot_ const &first_knot = GetFront(), &last_knot = GetBack();
  if (IsLess(parametric_coordinate, first_knot, tolerance) || IsGreater(parametric_coordinate, last_knot, tolerance))
      throw OutOfRange("The parametric coordinate " + to_string(parametric_coordinate.Get()) + " is outside of the "
                       "knot vector's scope [" + to_string(first_knot.Get()) + "," + to_string(last_knot.Get()) + "].");
}

void KnotVector::ThrowIfTooSmallOrNotNonDecreasing(Tolerance const &tolerance) const {
  int const &number_of_knots = knots_.size();
  if (number_of_knots < 2) throw DomainError("The knot vector has to contain at least 2 knots but only contains " +
                                             to_string(number_of_knots) + ".");
  Index::ForEach(1, number_of_knots, [&] (Index const &knot) {
      Index::Type_ const &index = knot.Get();
      Knot_ const &current_knot = knots_[index], &previous_knot = knots_[index - 1];
      if (IsLess(current_knot, previous_knot, tolerance))  // NOLINT(readability/braces)
          throw DomainError("The knot vector has to be a non-decreasing sequence of real numbers but the knot " +
                            to_string(current_knot.Get()) + " at index " + to_string(index) + " is less than the "
                            "previous knot " + to_string(previous_knot.Get()) + "."); });
}
#endif

}  // namespace splinelib::sources::parameter_spaces
