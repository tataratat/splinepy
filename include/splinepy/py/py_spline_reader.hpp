#pragma once

#include <string>

// pybind11
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

// SplineLib
#include <BSplineLib/InputOutput/iges.hpp>
#include <BSplineLib/Splines/b_spline.hpp>
#include <BSplineLib/Utilities/named_type.hpp>

namespace splinepy::py {

namespace py = pybind11;

using namespace bsplinelib;

/// @brief BSpline parser class
/// @tparam para_dim Dimension of parametric space
/// @tparam dim Dimension of physical space
template<int para_dim, int dim>
class BSplineParser {
public:
  using BSpline = typename splines::BSpline<para_dim, dim>;
  using OutputInformation = typename BSpline::OutputInformation_;
  using OutputParameterSpace =
      typename std::tuple_element_t<0, OutputInformation>;
  using OutputVectorSpace = typename std::tuple_element_t<1, OutputInformation>;
  using OutputKnotVectors =
      typename std::tuple_element_t<0, OutputParameterSpace>;
  using OutputCoordinates = typename std::tuple_element_t<0, OutputVectorSpace>;
  using OutputDegrees = typename std::tuple_element_t<1, OutputParameterSpace>;

  using SplineEntry = input_output::SplineEntry;

  /// @brief Convert BSpline to Python dictionary
  /// @param bspline BSpline object
  static py::dict BSplineToDict(SplineEntry const& bspline);
};

/// @brief BSpline parser class
/// @tparam para_dim Dimension of parametric space
/// @tparam dim Dimension of physical space
template<int para_dim, int dim>
class NurbsParser {
public:
  using Nurbs = typename splines::Nurbs<para_dim, dim>;
  using OutputInformation = typename Nurbs::OutputInformation_;
  using OutputParameterSpace =
      typename std::tuple_element_t<0, OutputInformation>;
  using OutputWeightedVectorSpace =
      typename std::tuple_element_t<1, OutputInformation>;
  using OutputKnotVectors =
      typename std::tuple_element_t<0, OutputParameterSpace>;
  using OutputCoordinates =
      typename std::tuple_element_t<0, OutputWeightedVectorSpace>;
  using OutputWeights =
      typename std::tuple_element_t<1, OutputWeightedVectorSpace>;
  using OutputDegrees = typename std::tuple_element_t<1, OutputParameterSpace>;

  using SplineEntry = input_output::SplineEntry;

  /// @brief Convert Nurbs to Python dictionary
  /// @param nurbs Nurbs object
  static py::dict NurbsToDict(SplineEntry const& nurbs);
};

/// @brief Spline reader class
class SplineReader {
public:
  using Splines = input_output::Splines;
  using SplineEntry = input_output::SplineEntry;

  /// @brief Parse BSpline
  /// @param bspline BSpline object
  static py::dict ParseBSpline(SplineEntry const& bspline);

  /// @brief Parse Nurbs
  /// @param nurbs Nurbs object
  static py::dict ParseNurbs(SplineEntry const& nurbs);

  /// @brief Read spline in .iges form.
  /// @param fname File name
  static py::list ReadIges(std::string fname);
};

} // namespace splinepy::py
