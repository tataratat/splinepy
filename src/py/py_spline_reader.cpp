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
