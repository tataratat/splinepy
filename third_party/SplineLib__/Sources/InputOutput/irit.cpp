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

#include "Sources/InputOutput/irit.hpp"

#include <iterator>
#include <utility>

#include "Sources/Splines/spline.hpp"
#include "Sources/Utilities/error_handling.hpp"
#include "Sources/Utilities/std_container_operations.hpp"
#include "Sources/Utilities/string_operations.hpp"
#include "Sources/Utilities/system_operations.hpp"

namespace splinelib::sources::input_output::irit {

namespace {

template<int parametric_dimensionality>
using ParameterSpace = typename splines::Spline<parametric_dimensionality, parametric_dimensionality>::ParameterSpace_;
using StringVectorConstIterator = utilities::string_operations::StringVectorConstIterator;
using std::make_shared, std::move, std::to_string, utilities::string_operations::ConvertToNumber,
      utilities::string_operations::StartsWith, utilities::string_operations::TrimCharacter,
      utilities::system_operations::InputStream, utilities::system_operations::Open,
      utilities::system_operations::OutputStream;

void SkipAttributes(StringVectorConstIterator &entry);
template<int parametric_dimensionality>
SplineEntry CreateSpline(StringVectorConstIterator &entry);
template<int parametric_dimensionality, int dimensionality>
SplineEntry CreateSpline(SharedPointer<ParameterSpace<parametric_dimensionality>> parameter_space,
                         bool const &is_non_rational, StringVectorConstIterator &entry);

template<int parametric_dimensionality>
void WriteSpline(OutputStream &file, SplineEntry const &spline, String const &type, Precision const &precision);
template<int parametric_dimensionality, int dimensionality>
void WriteSpline(OutputStream &file, SplineEntry const &spline, Precision const &precision);
template<bool is_rational, typename SplineType>
void WriteSpline(OutputStream &file, SplineType const &spline, String const &point_type, Precision const &precision);

}  // namespace

Splines Read(String const &file_name) {
  Splines splines;
#ifndef NDEBUG
  try {
#endif
    InputStream file{Open<InputStream, utilities::system_operations::kModeIn>(file_name)};
    String line, content;
    while (getline(file, line)) content += line;
    StringVector const &entries = utilities::string_operations::SplitAtDelimiter(content, ' ');
    StringVectorConstIterator entry{entries.begin()};
    while (entry != std::prev(entries.end())) {
      String const &geometrical_type = TrimCharacter(*(entry++), '[');
      if (*entry == "BSPLINE") {
        SkipAttributes(++entry);
        if (geometrical_type == "CURVE") {
          splines.emplace_back(CreateSpline<1>(entry));
        } else if (geometrical_type == "SURFACE") {
          splines.emplace_back(CreateSpline<2>(entry));
        } else if (geometrical_type == "TRIVAR") {
          splines.emplace_back(CreateSpline<3>(entry));
        } else {
#ifndef NDEBUG
          throw RuntimeError("Unknown geometrical type (" + geometrical_type + ") encountered as currently only "
                             "entries of geometrical type CURVE, SURFACE, and TRIVAR BSPLINE can be read.");
#endif
        }
      }
    }
#ifndef NDEBUG
  } catch (RuntimeError const &exception) { Throw(exception, "splinelib::sources::input_output::irit::Read"); }
#endif
  return splines;
}

void Write(Splines const &splines, String const &file_name, Precision const &precision) {
#ifndef NDEBUG
  try {
#endif
    OutputStream file{Open<OutputStream, utilities::system_operations::kModeOut>(file_name)};
    file << "[OBJECT SPLINES";
    Index::ForEach(0, splines.size(), [&] (Index const &spline_index) {
        Index::Type_ const &spline_index_value = spline_index.Get();
        file << "\n  [OBJECT SPLINE_" << to_string(spline_index_value);
        SplineEntry const &spline = splines[spline_index_value];
        int const &parametric_dimensionality = spline->parametric_dimensionality_;
        switch (parametric_dimensionality) {
          case 1:
            WriteSpline<1>(file, spline, "CURVE", precision);
            break;
          case 2:
            WriteSpline<2>(file, spline, "SURFACE", precision);
            break;
          case 3:
            WriteSpline<3>(file, spline, "TRIVAR", precision);
            break;
          default:
#ifndef NDEBUG
            throw RuntimeError("The parametric dimensionality (" + to_string(parametric_dimensionality) + ") of the "
                      "spline number " + to_string(spline_index_value) + " must be greater than 0 and less than 4.");
#endif
            break;
        }
        file << "\n  ]\n"; });
    file << "]";
#ifndef NDEBUG
  } catch (RuntimeError const &exception) { Throw(exception, "splinelib::sources::input_output::irit::Write"); }
#endif
}

namespace {

void SkipAttributes(StringVectorConstIterator &entry) {
  while (StartsWith(*entry, "[")) while (!utilities::string_operations::EndsWith(*(entry++), "]")) {}
}

template<int parametric_dimensionality>
SplineEntry CreateSpline(StringVectorConstIterator &entry) {
  using DimensionType = Dimension::Type_;
  using ParameterSpace = ParameterSpace<parametric_dimensionality>;
  using KnotVectors = typename ParameterSpace::KnotVectors_;
  using KnotVector = typename KnotVectors::value_type::element_type;
  using Knot = typename KnotVector::Knot_;

  typename ParameterSpace::NumberOfBasisFunctions_ number_of_coordinates;
  Dimension::ForEach(0, parametric_dimensionality, [&] (Dimension const &dimension) {
      number_of_coordinates[dimension.Get()] = ConvertToNumber<Length>(*(entry++)); });
  typename ParameterSpace::Degrees_ degrees;
  Dimension::ForEach(0, parametric_dimensionality, [&] (Dimension const &dimension) {
      degrees[dimension.Get()] = --ConvertToNumber<Degree>(*(entry++)); });
  String const &point_type = *(entry++);
  KnotVectors knot_vectors;
  Dimension::ForEach(0, parametric_dimensionality, [&] (Dimension const &dimension) {
      DimensionType const &current_dimension = dimension.Get();
      SkipAttributes(++entry);
      int const &number_of_knots = number_of_coordinates[current_dimension].Get() + degrees[current_dimension].Get() +
                                   1;
      typename KnotVector::Knots_ knots{};
      knots.reserve(number_of_knots);
      knots.emplace_back(ConvertToNumber<Knot>(*(entry++)));
      Index::ForEach(1, number_of_knots - 1, [&] (Index const &) {
          knots.emplace_back(ConvertToNumber<Knot>(*(entry++))); });
      knots.emplace_back(ConvertToNumber<Knot>(TrimCharacter(*(entry++), ']')));
      knot_vectors[current_dimension] = make_shared<KnotVector>(knots); });
  SharedPointer<ParameterSpace> parameter_space{make_shared<ParameterSpace>(move(knot_vectors), move(degrees))};
  SplineEntry spline;
  DimensionType const &dimensionality = ConvertToNumber<DimensionType>(point_type.substr(1, 1));
  bool const &is_non_rational = StartsWith(point_type, "E");
  switch (dimensionality) {
    case 1:
      spline = CreateSpline<parametric_dimensionality, 1>(move(parameter_space), is_non_rational, entry);
      break;
    case 2:
      spline = CreateSpline<parametric_dimensionality, 2>(move(parameter_space), is_non_rational, entry);
      break;
    case 3:
      spline = CreateSpline<parametric_dimensionality, 3>(move(parameter_space), is_non_rational, entry);
      break;
    case 4:
      spline = CreateSpline<parametric_dimensionality, 4>(move(parameter_space), is_non_rational, entry);
      break;
    default:
#ifndef NDEBUG
      throw RuntimeError("The spline's dimensionality " + to_string(dimensionality) + " must be greater than 0 and "
                         "less than 10.");
#endif
      break;
  }
  return spline;
}

template<int parametric_dimensionality, int dimensionality>
SplineEntry CreateSpline(SharedPointer<ParameterSpace<parametric_dimensionality>> parameter_space,
                         bool const &is_non_rational, StringVectorConstIterator &entry) {
  int const &total_number_of_coordinates = parameter_space->GetTotalNumberOfBasisFunctions(),
            &maximum_dimension = (dimensionality - 1);
  if (is_non_rational) {
    using BSpline = BSpline<parametric_dimensionality, dimensionality>;
    using VectorSpace = typename BSpline::VectorSpace_;

    typename VectorSpace::Coordinates_ coordinates{};
    coordinates.reserve(total_number_of_coordinates);
    Index::ForEach(0, total_number_of_coordinates, [&] (Index const &) {
        typename VectorSpace::Coordinate_ scalar_coordinates;
        scalar_coordinates[0] = ConvertToNumber<Coordinate>(TrimCharacter(*(entry++), '['));
        Dimension::ForEach(1, maximum_dimension, [&] (Dimension const &dimension) {
            scalar_coordinates[dimension.Get()] = ConvertToNumber<Coordinate>(*(entry++)); });
        scalar_coordinates[maximum_dimension] = ConvertToNumber<Coordinate>(TrimCharacter(*(entry++), ']'));
        coordinates.emplace_back(scalar_coordinates); });
    return make_shared<BSpline>(move(parameter_space), make_shared<VectorSpace>(move(coordinates)));
  } else {
    using Nurbs = Nurbs<parametric_dimensionality, dimensionality>;
    using WeightedVectorSpace = typename Nurbs::WeightedVectorSpace_;

    typename WeightedVectorSpace::Weights_ weights{};
    weights.reserve(total_number_of_coordinates);
    typename WeightedVectorSpace::Coordinates_ coordinates{};
    coordinates.reserve(total_number_of_coordinates);
    Index::ForEach(0, total_number_of_coordinates, [&] (Index const &) {
        weights.emplace_back(ConvertToNumber<Weight>(TrimCharacter(*(entry++), '[')));
        typename WeightedVectorSpace::Coordinate_ scalar_coordinates;
        Dimension::ForEach(0, maximum_dimension, [&] (Dimension const &dimension) {
            scalar_coordinates[dimension.Get()] = ConvertToNumber<Coordinate>(*(entry++)); });
        scalar_coordinates[maximum_dimension] = ConvertToNumber<Coordinate>(TrimCharacter(*(entry++), ']'));
        coordinates.emplace_back(scalar_coordinates); });
    return make_shared<Nurbs>(move(parameter_space), make_shared<WeightedVectorSpace>(move(coordinates),
                                                                                      move(weights)));
  }
}

template<int parametric_dimensionality>
void WriteSpline(OutputStream &file, SplineEntry const &spline, String const &type, Precision const &precision) {
  file << "\n    [" << type << " BSPLINE";
  int const &dimensionality = spline->dimensionality_;
  switch (dimensionality) {
    case 1:
      WriteSpline<parametric_dimensionality, 1>(file, spline, precision);
      break;
    case 2:
      WriteSpline<parametric_dimensionality, 2>(file, spline, precision);
      break;
    case 3:
      WriteSpline<parametric_dimensionality, 3>(file, spline, precision);
      break;
    case 4:
      WriteSpline<parametric_dimensionality, 4>(file, spline, precision);
      break;
    default:
#ifndef NDEBUG
      throw RuntimeError("The spline's dimensionality (" + to_string(dimensionality) + ") must be greater than 0 and "
                         "less than 10.");
#endif
      break;
  }
}

template<int parametric_dimensionality, int dimensionality>
void WriteSpline(OutputStream &file, SplineEntry const &spline, Precision const &precision) {
  using std::static_pointer_cast;

  String const &dimensionality_string = to_string(dimensionality);
  spline->is_rational_ ? WriteSpline<true>(file, *static_pointer_cast<Nurbs<parametric_dimensionality,
      dimensionality>>(spline), "P" + dimensionality_string, precision) : WriteSpline<false>(file,
          *static_pointer_cast<BSpline<parametric_dimensionality, dimensionality>>(spline), "E" + dimensionality_string,
              precision);
}

template<bool is_rational, typename SplineType>
void WriteSpline(OutputStream &file, SplineType const &spline, String const &point_type, Precision const &precision) {
  using OutputInformation = typename SplineType::OutputInformation_;
  using std::tuple_element_t;
  using ParameterSpace = tuple_element_t<0, OutputInformation>;
  using VectorSpace = tuple_element_t<1, OutputInformation>;
  using std::get, utilities::string_operations::Append;

  OutputInformation const &spline_written = spline.Write(precision);
  ParameterSpace const &parameter_space = get<0>(spline_written);
  String number_of_coordinates_orders_and_point_type{};
  Append(number_of_coordinates_orders_and_point_type, " ", get<2>(parameter_space));
  for (typename tuple_element_t<1, ParameterSpace>::value_type const &degree : get<1>(parameter_space))
      Append(number_of_coordinates_orders_and_point_type, " ",
             to_string(ConvertToNumber<typename SplineType::ParameterSpace_::Degrees_::value_type::Type_>(degree) + 1));
  Append(number_of_coordinates_orders_and_point_type, " ", point_type);
  file << number_of_coordinates_orders_and_point_type;
  for (typename tuple_element_t<0, ParameterSpace>::value_type const &knot_vector : get<0>(parameter_space)) {
    String knot_vectors{"KV"};
    Append(knot_vectors, " ", knot_vector);
    file << "\n      [" << knot_vectors << "]";
  }
  VectorSpace const &vector_space = get<1>(spline_written);
  tuple_element_t<0, VectorSpace> const &coordinates = get<0>(vector_space);
  Index::ForEach(0, coordinates.size(), [&] (Index const &coordinate) {
      Index::Type_ const &coordinate_index = coordinate.Get();
      String coordinate_string{};
      if constexpr (is_rational) Append(coordinate_string, " ", get<1>(vector_space)[coordinate_index]);
      Append(coordinate_string, " ", coordinates[coordinate_index]);
      file << "\n      [" << coordinate_string.erase(0, 1) << "]"; });
  file << "\n    ]";
}

}  // namespace

}  // namespace splinelib::sources::input_output::irit
