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

template<int parametric_dimensionality, int dimensionality, bool is_rational>
auto CastToSpline(SplineEntry const &spline_entry) {
  using std::dynamic_pointer_cast;

  if constexpr (is_rational) {
    return dynamic_pointer_cast<Nurbs<parametric_dimensionality, dimensionality>>(spline_entry);
  } else {
    return dynamic_pointer_cast<BSpline<parametric_dimensionality, dimensionality>>(spline_entry);
  }
}

template<typename Coordinate>
String WriteCoordinate3d(Coordinate const &coordinate, String const &delimiter) {
  String written{coordinate[0]};
  Dimension::ForEach(1, 3, [&] (Dimension const &dimension) {
      Dimension::Type_ const &dimension_value = dimension.Get();
      utilities::string_operations::Append(written, delimiter, dimension_value < static_cast<Type>(coordinate.size()) ?
                                                                   coordinate[dimension_value] : "0"); });
  return written;
}
