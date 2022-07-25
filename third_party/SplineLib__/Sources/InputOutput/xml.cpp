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

#include "Sources/InputOutput/xml.hpp"

#include <algorithm>
#include <iterator>
#include <utility>

#include <pugixml.hpp>  // Cf. documentation at <https://pugixml.org/docs/manual.html>.
#include "Sources/Splines/spline.hpp"
#include "Sources/Utilities/error_handling.hpp"
#include "Sources/Utilities/std_container_operations.hpp"
#include "Sources/Utilities/string_operations.hpp"

namespace splinelib::sources::input_output::xml {

using File = pugi::xml_document;
using Node = pugi::xml_node;

namespace {

template<int parametric_dimensionality>
using ParameterSpace = typename splines::Spline<parametric_dimensionality, parametric_dimensionality>::ParameterSpace_;
using std::for_each, std::to_string, utilities::string_operations::ConvertToNumber,
      utilities::string_operations::ConvertToNumbers, utilities::string_operations::TrimCharacter;

template<int parametric_dimensionality>
SplineEntry CreateSpline(Node const &spline_entry);
template<int parametric_dimensionality, int dimensionality>
SplineEntry CreateSpline(SharedPointer<ParameterSpace<parametric_dimensionality>> parameter_space,
                         Node const &spline_entry);

template<int parametric_dimensionality>
void Write(Node &spline_node, SplineEntry const &spline, Precision const &precision);
template<int parametric_dimensionality, int dimensionality>
void Write(Node &spline_node, SplineEntry const &spline, String const &names, Precision const &precision);
template<bool is_rational, typename SplineType>
void Write(Node &spline_node, SplineType const &spline, String const &names, Precision const &precision);

}  // namespace

Splines Read(String const &file_name) {
  Splines splines;
#ifndef NDEBUG
  try {
#endif
    File file;
    pugi::xml_parse_result const &result = file.load_file(file_name.c_str());
#ifndef NDEBUG
    if (!result) throw RuntimeError("XML file " + file_name + " cannot be parsed.");
#endif
    Node const &spline_list = file.child("SplineList");
    for_each(spline_list.begin(), spline_list.end(), [&] (Node const &spline_entry) {
        int const &parametric_dimensionality = ConvertToNumber<int>(TrimCharacter(
                                                                        spline_entry.attribute("splDim").value(), ' '));
        SplineEntry spline;
        switch (parametric_dimensionality) {
          case 1:
            spline = CreateSpline<1>(spline_entry);
            break;
          case 2:
            spline = CreateSpline<2>(spline_entry);
            break;
          case 3:
            spline = CreateSpline<3>(spline_entry);
            break;
          case 4:
            spline = CreateSpline<4>(spline_entry);
            break;
          default:
#ifndef NDEBUG
            throw RuntimeError("The parametric dimensionality (" + to_string(parametric_dimensionality) + ") must be "
                               "greater than 0 and currently less than 5.");
#endif
            break;
        }
        splines.emplace_back(spline); });
#ifndef NDEBUG
  } catch (RuntimeError const &exception) { Throw(exception, "splinelib::sources::input_output::xml::Read"); }
#endif
  return splines;
}

void Write(Splines const &splines, String const &file_name, Precision const &precision) {
#ifndef NDEBUG
  try {
#endif
    File file;
    Node spline_list{file.append_child("SplineList")};
    spline_list.append_attribute("NumberOfSplines") = to_string(splines.size()).c_str();
    for_each(splines.begin(), splines.end(), [&] (SplineEntry const &spline) {
        Node spline_node{spline_list.append_child("SplineEntry")};
        int const &parametric_dimensionality = spline->parametric_dimensionality_;
        spline_node.append_attribute("splDim") = parametric_dimensionality;
        switch (parametric_dimensionality) {
          case 1:
            Write<1>(spline_node, spline, precision);
            break;
          case 2:
            Write<2>(spline_node, spline, precision);
            break;
          case 3:
            Write<3>(spline_node, spline, precision);
            break;
          case 4:
            Write<4>(spline_node, spline, precision);
            break;
          default:
    #ifndef NDEBUG
            throw RuntimeError("The spline's parametric dimensionality (" + to_string(parametric_dimensionality) + ") "
                               "must be greater than 0 and currently less than 5.");
    #endif
            break;
        } });
    file.save_file(file_name.c_str(), "  ", pugi::format_indent_attributes);
#ifndef NDEBUG
  } catch (RuntimeError const &exception) { Throw(exception, "splinelib::sources::input_output::xml::Write"); }
#endif
}

namespace {

using std::make_shared, std::move;

template<int parametric_dimensionality>
SplineEntry CreateSpline(Node const &spline_entry) {
  using ParameterSpace = ParameterSpace<parametric_dimensionality>;
  using Degrees = typename ParameterSpace::Degrees_;
  using KnotVectors = typename ParameterSpace::KnotVectors_;
  using KnotVector = typename KnotVectors::value_type::element_type;

  Node child{spline_entry.child("kntVecs").first_child()};
  KnotVectors knot_vectors;
  Dimension::ForEach(0, parametric_dimensionality, [&] (Dimension const &dimension) {
      knot_vectors[dimension.Get()] =
          make_shared<KnotVector>(ConvertToNumbers<typename KnotVector::Knot_>(child.first_child().value(), ' '));
      child = child.next_sibling(); });
  SharedPointer<ParameterSpace> parameter_space{make_shared<ParameterSpace>(move(knot_vectors),
      move(utilities::std_container_operations::TransformNamedTypes<Degrees>(
               ConvertToNumbers<typename Degrees::value_type>(spline_entry.child("deg").first_child().value(), ' '))))};
  SplineEntry spline;
  int const &dimensionality = ConvertToNumber<int>(TrimCharacter(spline_entry.attribute("spaceDim").value(), ' '));
  switch (dimensionality) {
    case 1:
      spline = CreateSpline<parametric_dimensionality, 1>(move(parameter_space), spline_entry);
      break;
    case 2:
      spline = CreateSpline<parametric_dimensionality, 2>(move(parameter_space), spline_entry);
      break;
    case 3:
      spline = CreateSpline<parametric_dimensionality, 3>(move(parameter_space), spline_entry);
      break;
    case 4:
      spline = CreateSpline<parametric_dimensionality, 4>(move(parameter_space), spline_entry);
      break;
    default:
#ifndef NDEBUG
      throw RuntimeError("The spline's dimensionality (" + to_string(dimensionality) + ") must be greater than 0 and "
                         "currently less than 10.");
#endif
      break;
  }
  return spline;
}

template<int parametric_dimensionality, int dimensionality>
SplineEntry CreateSpline(SharedPointer<ParameterSpace<parametric_dimensionality>> parameter_space,
                         Node const &spline_entry) {
  using BSpline = BSpline<parametric_dimensionality, dimensionality>;
  using Nurbs = Nurbs<parametric_dimensionality, dimensionality>;
  using VectorSpace = typename BSpline::VectorSpace_;
  using WeightedVectorSpace = typename Nurbs::WeightedVectorSpace_;
  using Coordinate = typename VectorSpace::Coordinate_;
  using ScalarCoordinate = typename Coordinate::value_type;
  using Coordinates = utilities::string_operations::Numbers<ScalarCoordinate>;

  StringVector const &names =
      utilities::string_operations::SplitAtDelimiter(spline_entry.child("cntrlPntVarNames").first_child().value(), ' ');
  utilities::string_operations::StringVectorConstIterator const &names_begin = names.begin(), &names_end = names.end(),
                                                                &x = std::find(names_begin, names.end(), "x");
  int const &start = (x != names_end ? static_cast<int>(std::distance(names_begin, x)) : 0);
  Coordinates const &coordinates_raw =
      ConvertToNumbers<ScalarCoordinate>(spline_entry.child("cntrlPntVars").first_child().value(), ' ');
  typename Coordinates::const_iterator scalar_coordinate_raw{coordinates_raw.begin() + start};
  int const &total_number_of_coordinates =
      ConvertToNumber<int>(TrimCharacter(spline_entry.attribute("numCntrlPnts").value(), ' '));
  typename VectorSpace::Coordinates_ coordinates{};
  coordinates.reserve(total_number_of_coordinates);
  Index::ForEach(0, total_number_of_coordinates, [&] (Index const &) {
      Coordinate scalar_coordinates;
      Dimension::ForEach(0, dimensionality, [&] (Dimension const &dimension) {
          scalar_coordinates[dimension.Get()] = ScalarCoordinate{*(scalar_coordinate_raw++)}; });
      coordinates.emplace_back(scalar_coordinates);
      scalar_coordinate_raw += (ConvertToNumber<int>(TrimCharacter(spline_entry.attribute("numOfCntrlPntVars").value(),
                                                                   ' ')) - dimensionality); });
  if (spline_entry.child("wght").empty()) {
    return make_shared<BSpline>(move(parameter_space), make_shared<VectorSpace>(move(coordinates)));
  } else {
    return make_shared<Nurbs>(move(parameter_space), make_shared<WeightedVectorSpace>(move(coordinates),
                                  move(ConvertToNumbers<typename WeightedVectorSpace::Weights_::value_type>(
                                           spline_entry.child("wght").first_child().value(), ' '))));
  }
}

template<int parametric_dimensionality>
void Write(Node &spline_node, SplineEntry const &spline, Precision const &precision) {
  int const &dimensionality = spline->dimensionality_;
  spline_node.append_attribute("spaceDim") = dimensionality;
  spline_node.append_attribute("numOfCntrlPntVars") = dimensionality;
  switch (dimensionality) {
    case 1:
      Write<parametric_dimensionality, 1>(spline_node, spline, "x", precision);
      break;
    case 2:
      Write<parametric_dimensionality, 2>(spline_node, spline, "x y", precision);
      break;
    case 3:
      Write<parametric_dimensionality, 3>(spline_node, spline, "x y z", precision);
      break;
    case 4:
      Write<parametric_dimensionality, 4>(spline_node, spline, "x y z t", precision);
      break;
    default:
#ifndef NDEBUG
      throw RuntimeError("The spline's dimensionality (" + to_string(dimensionality) + ") must be greater than 0 and "
                         "currently less than 5.");
#endif
      break;
  }
}

template<int parametric_dimensionality, int dimensionality>
void Write(Node &spline_node, SplineEntry const &spline, String const &names, Precision const &precision) {
  using std::static_pointer_cast;

  spline->is_rational_ ? Write<true>(spline_node, *static_pointer_cast<Nurbs<parametric_dimensionality,
      dimensionality>>(spline), names, precision) : Write<false>(spline_node,
          *static_pointer_cast<BSpline<parametric_dimensionality, dimensionality>>(spline), names, precision);
}

template<bool is_rational, typename SplineType>
void Write(Node &spline_node, SplineType const &spline, String const &names, Precision const &precision) {
  using OutputInformation = typename SplineType::OutputInformation_;
  using std::tuple_element_t;
  using ParameterSpace = tuple_element_t<0, OutputInformation>;
  using VectorSpace = tuple_element_t<1, OutputInformation>;
  using KnotVector = typename tuple_element_t<0, ParameterSpace>::value_type;
  using Coordinates = tuple_element_t<0, VectorSpace>;
  using std::get;

  OutputInformation const &spline_written = spline.Write(precision);
  VectorSpace const &vector_space = get<1>(spline_written);
  Coordinates const &coordinates = get<0>(vector_space);
  spline_node.append_attribute("numCntrlPnts") = coordinates.size();
  spline_node.append_child("cntrlPntVarNames").text() = ("\n      " + names + "\n    ").c_str();
  String coordinates_written;
  for_each(coordinates.begin(), coordinates.end(), [&] (typename Coordinates::value_type const &coordinate) {
      utilities::string_operations::Append(coordinates_written += "\n     ", " ", coordinate); });
  spline_node.append_child("cntrlPntVars").text() = (coordinates_written + "\n    ").c_str();
  if constexpr (is_rational) {
    String weights;
    for (typename tuple_element_t<1, VectorSpace>::value_type const &weight : get<1>(vector_space))
        weights += ("\n      " + weight);
    spline_node.append_child("wght").text() = (weights + "\n    ").c_str();
  }
  ParameterSpace const &parameter_space = get<0>(spline_written);
  String degrees;
  for (typename tuple_element_t<1, ParameterSpace>::value_type const &degree : get<1>(parameter_space))
      degrees += ("\n      " + degree);
  spline_node.append_child("deg").text() = (degrees + "\n    ").c_str();
  Node knot_vectors{spline_node.append_child("kntVecs")};
  for (KnotVector const &knot_vector : get<0>(parameter_space)) {
    String knot_vector_written{};
    for_each(knot_vector.begin(), knot_vector.end(), [&] (typename KnotVector::value_type const &knot) {
        knot_vector_written += ("\n        " + knot); });
    knot_vectors.append_child("kntVec").text() = (knot_vector_written + "\n      ").c_str();
  }
}

}  // namespace

}  // namespace splinelib::sources::input_output::xml
