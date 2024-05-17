/*
MIT License

Copyright (c) 2021 Jaewook Lee

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <memory>

#include <BSplineLib/InputOutput/iges.hpp>
#include <BSplineLib/InputOutput/operations.hpp>
#include <BSplineLib/Splines/spline_item.hpp>

#include "splinepy/py/py_spline.hpp"
#include "splinepy/py/py_spline_reader.hpp"
#include "splinepy/splines/bspline.hpp"
#include "splinepy/splines/nurbs.hpp"
#include "splinepy/utils/print.hpp"

namespace splinepy::py {

template<int para_dim>
std::shared_ptr<PySpline>
ConvertToPySpline(std::shared_ptr<bsplinelib::splines::SplineItem>& spline) {
  // rational -> NURBS
  if (spline->is_rational_) {
    auto bsplinelib_nurbs =
        bsplinelib::input_output::operations::CastToSpline<para_dim, true>(
            spline);
    auto splinepy_nurbs =
        std::make_shared<splinepy::splines::Nurbs<para_dim>>();
    // with proper default init, we can just take members as they are shared
    splinepy_nurbs->ShareMembers(bsplinelib_nurbs);

    return std::make_shared<PySpline>(splinepy_nurbs);
  }

  auto bsplinelib_bspline =
      bsplinelib::input_output::operations::CastToSpline<para_dim, false>(
          spline);
  auto splinepy_bspline =
      std::make_shared<splinepy::splines::BSpline<para_dim>>();
  splinepy_bspline->ShareMembers(bsplinelib_bspline);

  // else -> BSpline
  return std::make_shared<PySpline>(splinepy_bspline);
}

std::shared_ptr<PySpline>
ToPySpline(std::shared_ptr<bsplinelib::splines::SplineItem>& spline) {
  switch (spline->parametric_dimensionality_) {
  case 1:
    return ConvertToPySpline<1>(spline);
  case 2:
    return ConvertToPySpline<2>(spline);
  default:
    splinepy::utils::PrintAndThrowError(
        "invalid para_dim (",
        spline->parametric_dimensionality_,
        ") detected while loading spline. It should be less then 3.");
    break;
  }
  return std::shared_ptr<PySpline>{};
}

py::list ReadIges(const std::string fname) {
  // read
  py::list splines;

  // load
  auto c_splines = bsplinelib::input_output::iges::Read(fname);

  // convert
  for (std::shared_ptr<bsplinelib::splines::SplineItem>& spline : c_splines) {
    splines.append(ToPySpline(spline)->ToDerived());
  }

  return splines;
}

///  @brief Adds spline reader. Keys are
/// ["knot_vectors", "control_points", "degrees"] (+ ["weights"])
/// @param m Python module
void init_spline_reader(py::module_& m) {
  m.def("read_iges", &ReadIges, py::arg("fname"));
}

} // namespace splinepy::py
