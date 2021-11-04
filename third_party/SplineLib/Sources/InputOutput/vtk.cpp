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

#include "Sources/InputOutput/vtk.hpp"

#include <algorithm>

#include "Sources/Splines/b_spline.hpp"
#include "Sources/Splines/nurbs.hpp"
#include "Sources/Utilities/error_handling.hpp"
#include "Sources/Utilities/index.hpp"
#include "Sources/Utilities/std_container_operations.hpp"
#include "Sources/Utilities/string_operations.hpp"
#include "Sources/Utilities/system_operations.hpp"

namespace splinelib::sources::input_output::vtk {

namespace {

template<int parametric_dimensionality>
using Index = utilities::Index<parametric_dimensionality>;
using NumberOfParametricCoordinates = NumbersOfParametricCoordinates::value_type;
using ScalarIndex = splinelib::Index;
using std::to_string, utilities::string_operations::Append;

template<int parametric_dimensionality>
void Sample(SplineEntry const &spline, NumberOfParametricCoordinates const &number_of_parametric_coordinates,
            int const &number_of_vertices_per_cell, String const &cell_type, int &total_number_of_points,
            String &points, int &total_number_of_cells, String &cells, int &cell_list_size, String &cell_types,
            Tolerance const &tolerance, Precision const &precision);
template<int parametric_dimensionality, int dimensionality>
void Sample(SplineEntry const &spline, NumberOfParametricCoordinates const &number_of_parametric_coordinates,
            int const &number_of_vertices_per_cell, String const &cell_type,  int &total_number_of_points,
            String &points, int &total_number_of_cells, String &cells, int &cell_list_size, String &cell_types,
            Tolerance const &tolerance, Precision const &precision);
template<int parametric_dimensionality, int dimensionality, typename SplineType>
void Sample(SplineType const &spline, NumberOfParametricCoordinates const &number_of_parametric_coordinates,
            int const &number_of_vertices_per_cell, String const &cell_type, int &total_number_of_points,
            String &points, int &total_number_of_cells, String &cells, int &cell_list_size, String &cell_types,
            Tolerance const &tolerance, Precision const &precision);

}  // namespace

void Sample(Splines const &splines, String const &file_name,
            NumbersOfParametricCoordinates const &numbers_of_parametric_coordinates, Tolerance const &tolerance,
            Precision const &precision) {
  using std::get, utilities::system_operations::OutputStream;

#ifndef NDEBUG
  try {
#endif
    int const &number_of_splines = splines.size(),
              &numbers_of_parametric_coordinates_size = numbers_of_parametric_coordinates.size();
#ifndef NDEBUG
    if (numbers_of_parametric_coordinates_size != number_of_splines)
        throw RuntimeError(to_string(number_of_splines) + " splines were given, but " +
            to_string(numbers_of_parametric_coordinates_size) + " numbers of parametric coordinates were provided.");
#endif
    OutputStream file{utilities::system_operations::Open<OutputStream,
                                                         utilities::system_operations::kModeOut>(file_name)};
    file << "# vtk DataFile Version 3.0\nSplines from Splinelib\nASCII\n\nDATASET UNSTRUCTURED_GRID";
    int total_number_of_points{}, total_number_of_cells{}, cell_list_size{};
    String points{}, cells{}, cell_types{};
    ScalarIndex::ForEach(0, number_of_splines, [&] (ScalarIndex const &spline_index) {
        ScalarIndex::Type_ const &spline_index_value = spline_index.Get();
        SplineEntry const &spline = splines[spline_index_value];
        int const &parametric_dimensionality = spline->parametric_dimensionality_;
        NumberOfParametricCoordinates const &number_of_parametric_coordinates =
            numbers_of_parametric_coordinates[spline_index_value];
#ifndef NDEBUG
        Message const spline_string{"For spline number " + to_string(spline_index_value) + ": "};

        int const &number_of_parametric_coordinates_size = number_of_parametric_coordinates.size();
        try {
          if (number_of_parametric_coordinates_size != parametric_dimensionality)
              throw RuntimeError("For splines of parametric dimensionality " +
                        to_string(parametric_dimensionality) + " numbers of parametric coordinates must be chosen, but "
                        "only " + to_string(number_of_parametric_coordinates_size) + " were provided.");
#endif
          switch (parametric_dimensionality) {
            case 1:
              Sample<1>(spline, number_of_parametric_coordinates, 2, "3", total_number_of_points, points,
                        total_number_of_cells, cells, cell_list_size, cell_types, tolerance, precision);
              break;
            case 2:
              Sample<2>(spline, number_of_parametric_coordinates, 4, "9", total_number_of_points, points,
                        total_number_of_cells, cells, cell_list_size, cell_types, tolerance, precision);
              break;
            case 3:
              Sample<3>(spline, number_of_parametric_coordinates, 8, "12", total_number_of_points, points,
                        total_number_of_cells, cells, cell_list_size, cell_types, tolerance, precision);
              break;
            default:
#ifndef NDEBUG
              throw RuntimeError("The spline's parametric dimensionality (" +
                        to_string(parametric_dimensionality) + ") must be larger than 0 and currently less than 4.");
#endif
              break;
          }
#ifndef NDEBUG
        } catch (RuntimeError const &exception) { throw RuntimeError(spline_string + exception.what()); }
#endif
    });
    file << "\nPOINTS " << total_number_of_points << " double" << points;
    file << "\n\nCELLS " << total_number_of_cells << " " << cell_list_size << cells;
    file << "\n\nCELL_TYPES " << total_number_of_cells << cell_types;
#ifndef NDEBUG
  } catch (RuntimeError const &exception) { Throw(exception, "splinelib::sources::input_output::vtk::Sample"); }
#endif
}

namespace {

template<int parametric_dimensionality>
void Sample(SplineEntry const &spline, NumberOfParametricCoordinates const &number_of_parametric_coordinates,
            int const &number_of_vertices_per_cell, String const &cell_type, int &total_number_of_points,
            String &points, int &total_number_of_cells, String &cells, int &cell_list_size, String &cell_types,
            Tolerance const &tolerance, Precision const &precision) {
  int const &dimensionality = spline->dimensionality_;
  switch (dimensionality) {
    case 1:
      Sample<parametric_dimensionality, 1>(spline, number_of_parametric_coordinates, number_of_vertices_per_cell,
                                           cell_type, total_number_of_points, points, total_number_of_cells, cells,
                                           cell_list_size, cell_types, tolerance, precision);
      break;
    case 2:
      Sample<parametric_dimensionality, 2>(spline, number_of_parametric_coordinates, number_of_vertices_per_cell,
                                           cell_type, total_number_of_points, points, total_number_of_cells, cells,
                                           cell_list_size, cell_types, tolerance, precision);
      break;
    case 3:
      Sample<parametric_dimensionality, 3>(spline, number_of_parametric_coordinates, number_of_vertices_per_cell,
                                           cell_type, total_number_of_points, points, total_number_of_cells, cells,
                                           cell_list_size, cell_types, tolerance, precision);
      break;
    default:
#ifndef NDEBUG
      throw RuntimeError("The spline's dimensionality (" + to_string(dimensionality) + ") must be larger than 0 and "
                         "less than 4.");
#endif
      break;
  }
}

template<int parametric_dimensionality, int dimensionality>
void Sample(SplineEntry const &spline, NumberOfParametricCoordinates const &number_of_parametric_coordinates,
            int const &number_of_vertices_per_cell, String const &cell_type, int &total_number_of_points,
            String &points, int &total_number_of_cells, String &cells, int &cell_list_size, String &cell_types,
            Tolerance const &tolerance, Precision const &precision) {
  using std::static_pointer_cast;

  if (spline->is_rational_) {
    Sample<parametric_dimensionality, dimensionality>(*static_pointer_cast<splines::Nurbs<parametric_dimensionality,
        dimensionality>>(spline), number_of_parametric_coordinates, number_of_vertices_per_cell, cell_type,
            total_number_of_points, points, total_number_of_cells, cells, cell_list_size, cell_types, tolerance,
                precision);
  } else {
    Sample<parametric_dimensionality, dimensionality>(*static_pointer_cast<splines::BSpline<parametric_dimensionality,
        dimensionality>>(spline), number_of_parametric_coordinates, number_of_vertices_per_cell, cell_type,
            total_number_of_points, points, total_number_of_cells, cells, cell_list_size, cell_types, tolerance,
                precision);
  }
}

template<int parametric_dimensionality, int dimensionality, typename SplineType>
void Sample(SplineType const &spline, NumberOfParametricCoordinates const &number_of_parametric_coordinates,
            int const &number_of_vertices_per_cell, String const &cell_type, int &total_number_of_points,
            String &points, int &total_number_of_cells, String &cells, int &cell_list_size, String &cell_types,
            Tolerance const &tolerance, Precision const &precision) {
  using Coordinates = typename SplineType::Base_::Coordinates_;
  using Index = Index<parametric_dimensionality>;
  using IndexLength = typename Index::Length_;
  using Vertices = Vector<ScalarIndex>;
  using std::for_each, utilities::string_operations::Write;

  constexpr Dimension const kDimension0{}, kDimension1{1};

  IndexLength number_of_vertices;
  std::copy(number_of_parametric_coordinates.begin(), number_of_parametric_coordinates.end(),
            number_of_vertices.begin());
  Coordinates const &current_points = spline.Sample(number_of_vertices, tolerance);
  for_each(current_points.begin(), current_points.end(), [&] (typename Coordinates::value_type const &point) {
    Append(points, "\n", operations::WriteCoordinate3d(Write<StringArray<dimensionality>>(point, precision), " "));
  });
  IndexLength number_of_cells;
  std::transform(number_of_parametric_coordinates.begin(), number_of_parametric_coordinates.end(),
      number_of_cells.begin(), [] (typename IndexLength::value_type sampling_value) { return --sampling_value; });
  String cells_string{};
  Index cell{Index::First(number_of_cells)};
  int const &current_total_number_of_cells = cell.GetTotalNumberOfIndices();
  cell_list_size += (current_total_number_of_cells * (number_of_vertices_per_cell + 1));
  Append(cell_types, "\n", StringVector(current_total_number_of_cells, cell_type));
  for (; cell != Index::Behind(number_of_cells); ++cell) {
    Index vertex{number_of_vertices, cell.GetIndex()};
    Vertices vertices;
    vertices.reserve(number_of_vertices_per_cell);
    vertices.emplace_back(vertex.GetIndex1d());
    vertices.emplace_back(vertex.Increment(kDimension0).GetIndex1d());
    if constexpr (parametric_dimensionality >= 2) {
      vertices.emplace_back(vertex.Increment(kDimension1).GetIndex1d());
      vertices.emplace_back(vertex.Decrement(kDimension0).GetIndex1d());
    }
    if constexpr (parametric_dimensionality == 3) {
      vertices.emplace_back(vertex.Decrement(kDimension1).Increment(Dimension{2}).GetIndex1d());
      vertices.emplace_back(vertex.Increment(kDimension0).GetIndex1d());
      vertices.emplace_back(vertex.Increment(kDimension1).GetIndex1d());
      vertices.emplace_back(vertex.Decrement(kDimension0).GetIndex1d());
    }
    Append(cells_string, "\n", to_string(number_of_vertices_per_cell));
    for_each(vertices.begin(), vertices.end(), [&] (typename Vertices::value_type const &vertex_index) {
        Append(cells_string, " ", to_string(total_number_of_points + vertex_index.Get())); });
  }
  total_number_of_points += current_points.size();
  total_number_of_cells += current_total_number_of_cells;
  Append(cells, "", cells_string);
}

}  // namespace

}  // namespace splinelib::sources::input_output::vtk
