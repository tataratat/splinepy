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

#include "Sources/InputOutput/iges.hpp"

#include <algorithm>
#include <limits>
#include <utility>

#include "Sources/Utilities/error_handling.hpp"
#include "Sources/Utilities/std_container_operations.hpp"
#include "Sources/Utilities/string_operations.hpp"
#include "Sources/Utilities/system_operations.hpp"

namespace splinelib::sources::input_output::iges {

namespace {

template<int parametric_dimensionality, int dimensionality>
using BSpline = BSpline<parametric_dimensionality, dimensionality>;
using utilities::string_operations::Numbers;
template<int parametric_dimensionality, int dimensionality>
using Nurbs = Nurbs<parametric_dimensionality, dimensionality>;
using SplineDataDouble = Numbers<double>;
using SplineDataInt = Numbers<int>;
using std::for_each, std::to_string, utilities::string_operations::Append,
      utilities::string_operations::ConvertToNumber, utilities::system_operations::InputStream,
      utilities::system_operations::Open, utilities::system_operations::OutputStream;

template<int parametric_dimensionality>
SplineEntry CreateSpline(SplineDataInt const &spline_data_int, SplineDataDouble const &spline_data_double,
                         int const &knot_vector_start);

int WriteSection(OutputStream &file, String const &section_content, String const &section_identifier);
template<int parametric_dimensionality>
Coordinate WriteSpline(SplineEntry const &spline, String const &delimiter, String &parameter_data_section_contribution,
                       Precision const &precision, Tolerance const &tolerance);
template<int parametric_dimensionality, int dimensionality>
Coordinate WriteSpline(SplineEntry const &spline, String const &delimiter, String &parameter_data_section_contribution,
                       Precision const &precision, Tolerance const &tolerance);
template<int parametric_dimensionality, int dimensionality, bool is_rational, typename SplineType>
void WriteSpline(SplineType const &spline, String const &delimiter, String &parameter_data_section_contribution,
                 Precision const &precision);
String ConvertToHollerith(String const &raw_string);
String GetLocalTime();
String MakeBlock(String const &raw_string, int const &block_width, bool const &align_right);

}  // namespace

Splines Read(String const &file_name) {
  using SplineSection = std::pair<Index, int>;
  using utilities::std_container_operations::GetValue, utilities::string_operations::ConvertToNumbers,
        utilities::string_operations::TrimCharacter;

  InputStream file;
#ifndef NDEBUG
  try {
#endif
    file = Open<InputStream, utilities::system_operations::kModeIn>(file_name);
#ifndef NDEBUG
  } catch (RuntimeError const &exception) {
    Throw(exception, "splinelib::sources::input_output::iges::Read");
  }
#endif
  String line;
  Vector<SplineSection> splines_read;
  // TODO(all): use all information from global section (e.g., delimiters, representatitivity, scaling & units)?
  do {} while (getline(file, line) && (GetValue(line, Index{72}) != 'D'));
  do {
    int const &entry_type = ConvertToNumber<int>(line.substr(5, 3));
    if ((entry_type == 126) || (entry_type == 128)) {
      Index const start{--ConvertToNumber<Index>(TrimCharacter(line.substr(8, 8), ' '))};
      getline(file, line);
      splines_read.emplace_back(start, ConvertToNumber<int>(TrimCharacter(line.substr(24, 8), ' ')));
    }
  } while (getline(file, line) && (GetValue(line, Index{72}) != 'P'));
  Splines splines;
  Index line_number{};
  for_each(splines_read.begin(), splines_read.end(), [&] (SplineSection const &spline) {
      Index const &spline_section_start = spline.first;
      for (; line_number < spline_section_start; getline(file, line), ++line_number) {}
      SplineDataInt const &spline_data_int = ConvertToNumbers<SplineDataInt::value_type>(
                                                 TrimCharacter(line.substr(0, 22), ','), ',');
      String spline_data_string{};
      for (; line_number < (spline_section_start + Index{spline.second}); getline(file, line), ++line_number)
          spline_data_string.append(TrimCharacter(line.substr(0, 64), ' '));
      SplineDataDouble const &spline_data_double = ConvertToNumbers<SplineDataDouble::value_type>(
                                                       TrimCharacter(spline_data_string, ';'),
                                                                            ',');
      if (*spline_data_int.begin() == 126) {
        splines.push_back(CreateSpline<1>(spline_data_int, spline_data_double, 7));
      } else {
        splines.push_back(CreateSpline<2>(spline_data_int, spline_data_double, 10));
      } });
  return splines;
}

void Write(Splines const &splines, String const &file_name, Precision const &precision, Tolerance const &tolerance) {
  using std::numeric_limits;

#ifndef NDEBUG
  try {
#endif
    OutputStream file{Open<OutputStream, utilities::system_operations::kModeOut>(file_name)};
    int const &start = WriteSection(file, MakeBlock("IGES file generated by SplineLib", 72, false), "S");
    Coordinate maximum_coordinate{};
    int directory_entry{1}, parameter_data{1};
    String entity_type, entity_name, directory_entry_section{}, parameter_data_section{};
    for_each(splines.begin(), splines.end(), [&] (SplineEntry const &spline) {
        int const &parametric_dimensionality = spline->parametric_dimensionality_;
        Coordinate maximum_coordinate_of_spline;
        String parameter_data_section_contribution{};
        switch (parametric_dimensionality) {
          case 1:
            parameter_data_section_contribution = ("126");
            maximum_coordinate_of_spline = WriteSpline<1>(spline, ",", parameter_data_section_contribution, precision,
                                                          tolerance);
            entity_type = MakeBlock("126", 8, true), entity_name = MakeBlock("BSp Crv", 8, true);
            break;
          case 2:
            parameter_data_section_contribution = ("128");
            maximum_coordinate_of_spline = WriteSpline<2>(spline, ",", parameter_data_section_contribution, precision,
                                                          tolerance);
            entity_type = MakeBlock("128", 8, true), entity_name = MakeBlock("BSp Surf", 8, true);
            break;
          default:
  #ifndef NDEBUG
            throw RuntimeError("The spline's parametric dimensionality (" + to_string(parametric_dimensionality) + ") "
                               "must be greater than 0 and less than 3.");
  #endif
            break;
        }
        if (maximum_coordinate <= maximum_coordinate_of_spline) maximum_coordinate = maximum_coordinate_of_spline;
        parameter_data_section_contribution += ";";
        int const &parameter_data_section_contribution_size = parameter_data_section_contribution.size();
        int parameter_data_section_contribution_number_of_lines{}, position{};
        while (position < parameter_data_section_contribution_size) {
          ++parameter_data_section_contribution_number_of_lines;
          int number_of_characters{};
          if ((position + 64) < parameter_data_section_contribution_size) {
            int character{};
            do {
              if (parameter_data_section_contribution[position + character++] == ',') number_of_characters = character;
            } while (character < 64);
          } else {
            number_of_characters = parameter_data_section_contribution_size - position;
          }
          parameter_data_section += (MakeBlock(parameter_data_section_contribution.substr(position,
              number_of_characters), 64, false) + MakeBlock(to_string(directory_entry), 8, true));
          position += number_of_characters;
        }
        String const &zero = MakeBlock("0", 8, true), &reserved = MakeBlock("", 8, false);
        Append(directory_entry_section, "", StringArray<18>{entity_type, MakeBlock(to_string(parameter_data), 8, true),
                   zero, zero, zero, zero, zero, zero, "00000001", entity_type, zero, zero, MakeBlock(to_string(
                       parameter_data_section_contribution_number_of_lines), 8, true), zero, reserved, reserved,
                           entity_name, zero});
        directory_entry += 2;
        parameter_data += parameter_data_section_contribution_number_of_lines; });
    String const &native_system_id = ConvertToHollerith("SplineLib 0.1.1"), &time = GetLocalTime();
    String global_section{ConvertToHollerith(",")};
    Append(global_section, ",", StringArray<25>{ConvertToHollerith(";"), ConvertToHollerith("product ID from sender"),
        ConvertToHollerith(file_name), native_system_id, native_system_id, to_string(numeric_limits<int>::digits10),
            to_string(numeric_limits<float>::max_exponent10), to_string(numeric_limits<float>::digits10),
                to_string(numeric_limits<double>::max_exponent10), to_string(numeric_limits<double>::digits10),
                    ConvertToHollerith("product ID for receiver"), "1.0", "6", ConvertToHollerith("M"), "1",
                        "0.0024384", time, "0.00000254", utilities::string_operations::Write(maximum_coordinate,
                            precision), ConvertToHollerith("Author"), ConvertToHollerith("Organization"), "11", "0",
                                time, ConvertToHollerith("NULL")});
    global_section += ";";
    int const &global = WriteSection(file, global_section, "G"), &directory = WriteSection(file,
        directory_entry_section, "D"), &parameter = WriteSection(file, parameter_data_section, "P");
    WriteSection(file, "S" + MakeBlock(to_string(start), 7, true) + "G" + MakeBlock(to_string(global), 7, true) +
                 "D" + MakeBlock(to_string(directory), 7, true) + "P" + MakeBlock(to_string(parameter), 7,
                     true) + MakeBlock("", 40, false), "T");
#ifndef NDEBUG
  } catch (RuntimeError const &exception) { Throw(exception, "splinelib::sources::input_output::iges::Write"); }
#endif
}

namespace {

template<int parametric_dimensionality>
SplineEntry CreateSpline(SplineDataInt const &spline_data_int, SplineDataDouble const &spline_data_double,
                         int const &knot_vector_start) {
  using BSpline = BSpline<parametric_dimensionality, 3>;
  using Nurbs = Nurbs<parametric_dimensionality, 3>;
  using ParameterSpace = typename BSpline::ParameterSpace_;
  using VectorSpace = typename BSpline::VectorSpace_;
  using WeightedVectorSpace = typename Nurbs::WeightedVectorSpace_;
  using KnotVectors = typename ParameterSpace::KnotVectors_;
  using KnotVector = typename KnotVectors::value_type::element_type;
  using std::make_shared, std::move;

  SplineDataInt::const_iterator spline_datum_int{spline_data_int.begin() + 1};
  typename ParameterSpace::NumberOfBasisFunctions_ number_of_coordinates;
  Dimension::ForEach(0, parametric_dimensionality, [&] (Dimension const &dimension) {
      number_of_coordinates[dimension.Get()] = Length{*(spline_datum_int++) + 1}; });
  typename ParameterSpace::Degrees_ degrees;
  Dimension::ForEach(0, parametric_dimensionality, [&] (Dimension const &dimension) {
      degrees[dimension.Get()] = Degree{*(spline_datum_int++)}; });
  SplineDataDouble::const_iterator spline_datum_double{spline_data_double.begin() + knot_vector_start};
  KnotVectors knot_vectors;
  Dimension::ForEach(0, parametric_dimensionality, [&] (Dimension const &dimension) {
      Dimension::Type_ const &current_dimension = dimension.Get();
      typename ParameterSpace::Knots_ knots{};
      Index::ForEach(0, number_of_coordinates[current_dimension].Get() + degrees[current_dimension].Get() + 1,
                     [&] (Index const &) { knots.emplace_back(*(spline_datum_double++)); });
      knot_vectors[current_dimension] = make_shared<KnotVector>(knots); });
  SharedPointer<ParameterSpace> parameter_space{make_shared<ParameterSpace>(move(knot_vectors), move(degrees))};
  int const &total_number_of_coordinates = parameter_space->GetTotalNumberOfBasisFunctions();
  typename WeightedVectorSpace::Weights_ weights{};
  weights.reserve(total_number_of_coordinates);
  Index::ForEach(0, total_number_of_coordinates,
                 [&] (Index const &) { weights.emplace_back(*(spline_datum_double++)); });
  typename VectorSpace::Coordinates_ coordinates{};
  coordinates.reserve(total_number_of_coordinates);
  Index::ForEach(0, total_number_of_coordinates, [&] (Index const &) {
      coordinates.push_back({Coordinate{*(spline_datum_double++)}, Coordinate{*(spline_datum_double++)},
                             Coordinate{*(spline_datum_double++)}}); });
  if (*(spline_datum_int + 2) == 1) {
    return make_shared<BSpline>(move(parameter_space), make_shared<VectorSpace>(move(coordinates)));
  } else {
    return make_shared<Nurbs>(move(parameter_space), make_shared<WeightedVectorSpace>(move(coordinates),
                                                                                      move(weights)));
  }
}

int WriteSection(OutputStream &file, String const &section_content, String const &section_identifier) {
  constexpr int const kBlockWidth{72};

  int line{};
  for (; line <= (static_cast<int>(section_content.size() - 1) / kBlockWidth); ) {
    file << MakeBlock(section_content.substr(line * kBlockWidth, kBlockWidth), kBlockWidth, false);
    file << (section_identifier + MakeBlock(to_string(++line), 7, true)) << std::endl;
  }
  return line;
}

template<int parametric_dimensionality>
Coordinate WriteSpline(SplineEntry const &spline, String const &delimiter, String &parameter_data_section_contribution,
                       Precision const &precision, Tolerance const &tolerance) {
  int const &dimensionality = spline->dimensionality_;
  switch (dimensionality) {
    case 1:
      return WriteSpline<parametric_dimensionality, 1>(spline, delimiter, parameter_data_section_contribution,
                                                       precision, tolerance);
      break;
    case 2:
      return WriteSpline<parametric_dimensionality, 2>(spline, delimiter, parameter_data_section_contribution,
                                                       precision, tolerance);
      break;
    case 3:
      return WriteSpline<parametric_dimensionality, 3>(spline, delimiter, parameter_data_section_contribution,
                                                       precision, tolerance);
      break;
    default:
#ifndef NDEBUG
      throw RuntimeError("The spline's dimensionality (" + to_string(dimensionality) + ") must be greater than 0 and "
                         "less than 4.");
#endif
      return Coordinate{};
      break;
  }
}

template<int parametric_dimensionality, int dimensionality>
Coordinate WriteSpline(SplineEntry const &spline, String const &delimiter, String &parameter_data_section_contribution,
                       Precision const &precision, Tolerance const &tolerance) {
  using std::static_pointer_cast;

  if (spline->is_rational_) {
    using Nurbs = Nurbs<parametric_dimensionality, dimensionality>;

    SharedPointer<Nurbs> const &nurbs = static_pointer_cast<Nurbs>(spline);
    WriteSpline<parametric_dimensionality, dimensionality, true>(*nurbs, delimiter, parameter_data_section_contribution,
                                                                 precision);
    return nurbs->ComputeUpperBoundForMaximumDistanceFromOrigin(tolerance);
  } else {
    using BSpline = BSpline<parametric_dimensionality, dimensionality>;

    SharedPointer<BSpline> const &b_spline = static_pointer_cast<BSpline>(spline);
    WriteSpline<parametric_dimensionality, dimensionality, false>(*b_spline, delimiter,
                                                                  parameter_data_section_contribution, precision);
    return b_spline->ComputeUpperBoundForMaximumDistanceFromOrigin(tolerance);
  }
}

template<int parametric_dimensionality, int dimensionality, bool is_rational, typename SplineType>
void WriteSpline(SplineType const &spline, String const &delimiter, String &parameter_data_section_contribution,
                 Precision const &precision) {
  using OutputInformation = typename SplineType::OutputInformation_;
  using std::tuple_element_t;
  using ParameterSpace = tuple_element_t<0, OutputInformation>;
  using VectorSpace = tuple_element_t<1, OutputInformation>;
  using KnotVectors = tuple_element_t<0, ParameterSpace>;
  using std::get;

  OutputInformation const &spline_written = spline.Write(precision);
  ParameterSpace const &parameter_space = get<0>(spline_written);
  for (typename tuple_element_t<2, ParameterSpace>::value_type const &number_of_basis_functions :
       get<2>(parameter_space))
      Append(parameter_data_section_contribution, delimiter, to_string(
                 ConvertToNumber<typename SplineType::ParameterSpace_::NumberOfBasisFunctions_::value_type::Type_>(
                     number_of_basis_functions) - 1));
  Append(parameter_data_section_contribution, delimiter, get<1>(parameter_space));
  Append(parameter_data_section_contribution, delimiter, StringArray<2>{"0", "0"});
  Append(parameter_data_section_contribution, delimiter, is_rational ? "0" : "1");
  Dimension::ForEach(0, parametric_dimensionality, [&] (Dimension const &) {
      Append(parameter_data_section_contribution, delimiter, "0"); });
  KnotVectors const &knot_vectors = get<0>(parameter_space);
  Append(parameter_data_section_contribution, delimiter, knot_vectors);
  VectorSpace const &vector_space = get<1>(spline_written);
  tuple_element_t<0, VectorSpace> const &coordinates = get<0>(vector_space);
  if constexpr (is_rational) {
    Append(parameter_data_section_contribution, delimiter, get<1>(vector_space));
  } else {
    Append(parameter_data_section_contribution, delimiter, StringVector(coordinates.size(), "1"));
  }
  for_each(coordinates.begin(), coordinates.end(),
      [&] (typename tuple_element_t<0, VectorSpace>::value_type const &coordinate) {
          Append(parameter_data_section_contribution, delimiter, operations::WriteCoordinate3d(coordinate, ",")); });
  Dimension::ForEach(0, parametric_dimensionality, [&] (Dimension const &dimension) {
      typename KnotVectors::value_type const &knot_vector = knot_vectors[dimension.Get()];
      Append(parameter_data_section_contribution, delimiter, knot_vector[0]);
      Append(parameter_data_section_contribution, delimiter, utilities::std_container_operations::GetBack(knot_vector));
  });
}

String ConvertToHollerith(String const &raw_string) {
  return (to_string(raw_string.size()) + "H" + raw_string);
}

String GetLocalTime() {
  utilities::system_operations::LocalTime const &local_time = utilities::system_operations::GetLocalTime();
  return ConvertToHollerith(to_string((local_time.tm_year + 1900) * 10000 + (local_time.tm_mon + 1) * 100 +
      local_time.tm_mday) + "." + to_string(local_time.tm_hour * 10000 + local_time.tm_min * 100 + local_time.tm_sec));
}

String MakeBlock(String const &raw_string, int const &block_width, bool const &align_right) {
  String const pad(block_width - static_cast<int>(raw_string.size()), ' ');
  return align_right ? pad + raw_string : raw_string + pad;
}

}  // namespace

}  // namespace splinelib::sources::input_output::iges
