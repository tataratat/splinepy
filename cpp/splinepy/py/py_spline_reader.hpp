#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <vector>

// pybind11
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

// SplineLib
#include <BSplineLib/InputOutput/iges.hpp>
#include <BSplineLib/Splines/b_spline.hpp>
#include <BSplineLib/Utilities/named_type.hpp>

#include <splinepy/utils/print.hpp>

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

  BSplineParser() {}

  /// @brief Convert BSpline to Python dictionary
  /// @param bspline BSpline object
  py::dict BSplineToDict(SplineEntry const& bspline) {
    int i, j;

    py::dict spline; //  to return

    std::shared_ptr<BSpline> bs = std::dynamic_pointer_cast<BSpline>(bspline);
    OutputInformation const& bs_info = bs->Write();

    //
    // Adapted from `bspline.hpp`
    //

    // Parameter space - knot vectors, degrees
    OutputParameterSpace const& parameter_space = std::get<0>(bs_info);
    OutputKnotVectors const& knot_vectors = std::get<0>(parameter_space);
    OutputDegrees const& degrees = std::get<1>(parameter_space);

    // Vector space - Coordinates(control points)
    OutputVectorSpace const& vector_space = std::get<1>(bs_info);
    OutputCoordinates const& coordinates = std::get<0>(vector_space);

    // Unpack - degrees
    auto p_degrees = py::array_t<int>(para_dim);
    py::buffer_info ds_buf = p_degrees.request();
    int* ds_buf_ptr = static_cast<int*>(ds_buf.ptr);

    i = 0;
    for (auto& degree : degrees) {
      ds_buf_ptr[i] =
          utilities::string_operations::ConvertToNumber<int>(degree);
      i++;
    }

    spline["degrees"] = p_degrees;

    // Unpack - knot vectors
    py::list p_knot_vectors;
    for (auto& knotvector : knot_vectors) {
      py::list p_kv;
      for (auto& knot : knotvector) {
        p_kv.append(
            utilities::string_operations::ConvertToNumber<double>(knot));
      }
      p_knot_vectors.append(p_kv);
    }

    spline["knot_vectors"] = p_knot_vectors;

    // Unpack - Coordinates (control points)
    auto p_control_points = py::array_t<double>(coordinates.size() * dim);
    py::buffer_info cps_buf = p_control_points.request();
    double* cps_buf_ptr = static_cast<double*>(cps_buf.ptr);

    i = 0;
    for (auto& coordinate : coordinates) {
      j = 0;
      for (auto& coord : coordinate) {
        cps_buf_ptr[i * dim + j] =
            utilities::string_operations::ConvertToNumber<double>(coord);
        j++;
      }
      i++;
    }

    p_control_points.resize({(int) coordinates.size(), dim});

    spline["control_points"] = p_control_points;

    return spline;
  }
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

  NurbsParser() {}

  /// @brief Convert Nurbs to Python dictionary
  /// @param nurbs Nurbs object
  py::dict NurbsToDict(SplineEntry const& nurbs) {
    int i, j;

    py::dict spline; //  to return

    // Read info
    std::shared_ptr<Nurbs> n = std::dynamic_pointer_cast<Nurbs>(nurbs);
    OutputInformation const& n_info = n->Write();

    //
    // Adapted from `nurbs.hpp`
    //

    // Parameter space - knot vectors, degrees
    OutputParameterSpace const& parameter_space = std::get<0>(n_info);
    OutputKnotVectors const& knot_vectors = std::get<0>(parameter_space);
    OutputDegrees const& degrees = std::get<1>(parameter_space);

    // Weighted vector space - Coordinates(control points), weights
    OutputWeightedVectorSpace const& weighted_vector_space =
        std::get<1>(n_info);
    OutputCoordinates const& coordinates = std::get<0>(weighted_vector_space);
    OutputWeights const& weights = std::get<1>(weighted_vector_space);

    // Unpack - Weights
    auto p_weights = py::array_t<double>(weights.size());
    py::buffer_info ws_buf = p_weights.request();
    double* ws_buf_ptr = static_cast<double*>(ws_buf.ptr);

    i = 0;
    for (auto& weight : weights) {
      ws_buf_ptr[i] =
          utilities::string_operations::ConvertToNumber<double>(weight);
      i++;
    }

    p_weights.resize({(int) weights.size(), 1}); // A tall vector

    spline["weights"] = p_weights;

    // Unpack - degrees
    auto p_degrees = py::array_t<int>(para_dim);
    py::buffer_info ds_buf = p_degrees.request();
    int* ds_buf_ptr = static_cast<int*>(ds_buf.ptr);

    i = 0;
    for (auto& degree : degrees) {
      ds_buf_ptr[i] =
          utilities::string_operations::ConvertToNumber<int>(degree);
      i++;
    }

    spline["degrees"] = p_degrees;

    // Unpack - knot vectors
    py::list p_knot_vectors;
    for (auto& knotvector : knot_vectors) {
      py::list p_kv;
      for (auto& knot : knotvector) {
        p_kv.append(
            utilities::string_operations::ConvertToNumber<double>(knot));
      }
      p_knot_vectors.append(p_kv);
    }

    spline["knot_vectors"] = p_knot_vectors;

    // Unpack - Coordinates (control points)
    auto p_control_points = py::array_t<double>(coordinates.size() * dim);
    py::buffer_info cps_buf = p_control_points.request();
    double* cps_buf_ptr = static_cast<double*>(cps_buf.ptr);

    i = 0;
    for (auto& coordinate : coordinates) {
      j = 0;
      for (auto& coord : coordinate) {
        cps_buf_ptr[i * dim + j] =
            utilities::string_operations::ConvertToNumber<double>(coord);
        j++;
      }
      i++;
    }

    p_control_points.resize({(int) coordinates.size(), dim});

    spline["control_points"] = p_control_points;

    return spline;
  }
};

/// @brief Spline reader class
class SplineReader {
public:
  using Splines = input_output::Splines;
  using SplineEntry = input_output::SplineEntry;

  /// @brief Spline object in C
  Splines c_splines;

  /// @brief Python-spline object
  py::list
      p_splines; // [dict(weights, degrees, knot_vectors, control_points), ...]

  /// @brief Read spline in .iges form.
  /// @param fname File name
  py::list ReadIges(std::string fname) {
    int i, j;

    c_splines = input_output::iges::Read(fname);
    Read(c_splines);

    return p_splines;
  }

  /// @brief Parse BSpline
  /// @param bspline BSpline object
  void ParseBSpline(SplineEntry const& bspline) {

    // Get paradim and dim
    int const& para_dim = bspline->parametric_dimensionality_;
    int const& dim = bspline->dimensionality_;

    // Find appropriate BSplines
    switch (para_dim) {
    case 1:
      switch (dim) {
      case 1: {
        auto bsparser = BSplineParser<1, 1>();
        p_splines.append(bsparser.BSplineToDict(bspline));
      } break;
      case 2: {
        auto bsparser = BSplineParser<1, 2>();
        p_splines.append(bsparser.BSplineToDict(bspline));
      } break;
      case 3: {
        auto bsparser = BSplineParser<1, 3>();
        p_splines.append(bsparser.BSplineToDict(bspline));
      } break;
      case 4: {
        auto bsparser = BSplineParser<1, 4>();
        p_splines.append(bsparser.BSplineToDict(bspline));
      } break;
      default:
        splinepy::utils::PrintAndThrowError(
            "ParseBSpline() does not support physical dimension of ",
            dim,
            ". Only dimensions 1 to 4 are supported.");
        break;
      }
      break;
    case 2:
      switch (dim) {
      case 1: {
        auto bsparser = BSplineParser<2, 1>();
        p_splines.append(bsparser.BSplineToDict(bspline));
      } break;
      case 2: {
        auto bsparser = BSplineParser<2, 2>();
        p_splines.append(bsparser.BSplineToDict(bspline));
      } break;
      case 3: {
        auto bsparser = BSplineParser<2, 3>();
        p_splines.append(bsparser.BSplineToDict(bspline));
      } break;
      case 4: {
        auto bsparser = BSplineParser<2, 4>();
        p_splines.append(bsparser.BSplineToDict(bspline));
      } break;
      default:
        splinepy::utils::PrintAndThrowError(
            "ParseBSpline() does not support physical dimension of ",
            dim,
            ". Only dimensions 1 to 4 are supported.");
        break;
      }
      break;
    case 3:
      switch (dim) {
      case 1: {
        auto bsparser = BSplineParser<3, 1>();
        p_splines.append(bsparser.BSplineToDict(bspline));
      } break;
      case 2: {
        auto bsparser = BSplineParser<3, 2>();
        p_splines.append(bsparser.BSplineToDict(bspline));
      } break;
      case 3: {
        auto bsparser = BSplineParser<3, 3>();
        p_splines.append(bsparser.BSplineToDict(bspline));
      } break;
      case 4: {
        auto bsparser = BSplineParser<3, 4>();
        p_splines.append(bsparser.BSplineToDict(bspline));
      } break;
      default:
        splinepy::utils::PrintAndThrowError(
            "ParseBSpline() does not support physical dimension of ",
            dim,
            ". Only dimensions 1 to 4 are supported.");
        break;
      }
      break;
    case 4:
      switch (dim) {
      case 1: {
        auto bsparser = BSplineParser<4, 1>();
        p_splines.append(bsparser.BSplineToDict(bspline));
      } break;
      case 2: {
        auto bsparser = BSplineParser<4, 2>();
        p_splines.append(bsparser.BSplineToDict(bspline));
      } break;
      case 3: {
        auto bsparser = BSplineParser<4, 3>();
        p_splines.append(bsparser.BSplineToDict(bspline));
      } break;
      case 4: {
        auto bsparser = BSplineParser<4, 4>();
        p_splines.append(bsparser.BSplineToDict(bspline));
      } break;
      default:
        splinepy::utils::PrintAndThrowError(
            "ParseBSpline() does not support physical dimension of ",
            dim,
            ". Only dimensions 1 to 4 are supported.");
        break;
      }
      break;
    default:
      splinepy::utils::PrintAndThrowError(
          "ParseBSpline() does not support parametric dimension of ",
          para_dim,
          ". Only dimensions 1 to 4 are supported.");
      break;
    }
  }

  /// @brief Parse Nurbs
  /// @param nurbs Nurbs object
  void ParseNurbs(SplineEntry const& nurbs) {

    // Get paradim and dim
    int const& para_dim = nurbs->parametric_dimensionality_;
    int const& dim = nurbs->dimensionality_;

    // Find appropriate Nurbs
    switch (para_dim) {
    case 1:
      switch (dim) {
      case 1: {
        auto nparser = NurbsParser<1, 1>();
        p_splines.append(nparser.NurbsToDict(nurbs));
      } break;
      case 2: {
        auto nparser = NurbsParser<1, 2>();
        p_splines.append(nparser.NurbsToDict(nurbs));
      } break;
      case 3: {
        auto nparser = NurbsParser<1, 3>();
        p_splines.append(nparser.NurbsToDict(nurbs));
      } break;
      case 4: {
        auto nparser = NurbsParser<1, 4>();
        p_splines.append(nparser.NurbsToDict(nurbs));
      } break;
      default:
        splinepy::utils::PrintAndThrowError(
            "ParseNurbs() does not support physical dimension of ",
            dim,
            ". Only dimensions 1 to 4 are supported.");
        break;
      }
      break;
    case 2:
      switch (dim) {
      case 1: {
        auto nparser = NurbsParser<2, 1>();
        p_splines.append(nparser.NurbsToDict(nurbs));
      } break;
      case 2: {
        auto nparser = NurbsParser<2, 2>();
        p_splines.append(nparser.NurbsToDict(nurbs));
      } break;
      case 3: {
        auto nparser = NurbsParser<2, 3>();
        p_splines.append(nparser.NurbsToDict(nurbs));
      } break;
      case 4: {
        auto nparser = NurbsParser<2, 4>();
        p_splines.append(nparser.NurbsToDict(nurbs));
      } break;
      default:
        splinepy::utils::PrintAndThrowError(
            "ParseNurbs() does not support physical dimension of ",
            dim,
            ". Only dimensions 1 to 4 are supported.");
        break;
      }
      break;
    case 3:
      switch (dim) {
      case 1: {
        auto nparser = NurbsParser<3, 1>();
        p_splines.append(nparser.NurbsToDict(nurbs));
      } break;
      case 2: {
        auto nparser = NurbsParser<3, 2>();
        p_splines.append(nparser.NurbsToDict(nurbs));
      } break;
      case 3: {
        auto nparser = NurbsParser<3, 3>();
        p_splines.append(nparser.NurbsToDict(nurbs));
      } break;
      case 4: {
        auto nparser = NurbsParser<3, 4>();
        p_splines.append(nparser.NurbsToDict(nurbs));
      } break;
      default:
        splinepy::utils::PrintAndThrowError(
            "ParseNurbs() does not support physical dimension of ",
            dim,
            ". Only dimensions 1 to 4 are supported.");
        break;
      }
      break;
    case 4:
      switch (dim) {
      case 1: {
        auto nparser = NurbsParser<4, 1>();
        p_splines.append(nparser.NurbsToDict(nurbs));
      } break;
      case 2: {
        auto nparser = NurbsParser<4, 2>();
        p_splines.append(nparser.NurbsToDict(nurbs));
      } break;
      case 3: {
        auto nparser = NurbsParser<4, 3>();
        p_splines.append(nparser.NurbsToDict(nurbs));
      } break;
      case 4: {
        auto nparser = NurbsParser<4, 4>();
        p_splines.append(nparser.NurbsToDict(nurbs));
      } break;
      default:
        splinepy::utils::PrintAndThrowError(
            "ParseNurbs() does not support physical dimension of ",
            dim,
            ". Only dimensions 1 to 4 are supported.");
        break;
      }
      break;
    default:
      splinepy::utils::PrintAndThrowError(
          "ParseNurbs() does not support parametric dimension of ",
          para_dim,
          ". Only dimensions 1 to 4 are supported.");
      break;
    }
  }

  /// @brief Reads spline
  /// @param splines
  void Read(Splines splines) {

    // Assign a new list
    //   - with `clear`, one object can't be reused: it alters all returned
    //   lists
    // Possible alternative is to return deepcopy of the list.
    //   - All the entries should be deepcopy-able
    p_splines = py::list();

    for (auto& spline : c_splines) {
      bool const& is_rational = spline->is_rational_;

      if (is_rational) {
        // Nurbs
        ParseNurbs(spline);
      } else {
        // BSpline
        ParseBSpline(spline);
      }
    }
  }
};

/* direct load calls */
/// Load iges
py::list ReadIges(std::string fname) {
  auto sr = SplineReader();
  return sr.ReadIges(fname);
}

///  @brief Adds spline reader. Keys are
/// ["knot_vectors", "control_points", "degrees"] (+ ["weights"])
/// @param m Python module
inline void add_spline_reader(py::module& m) {
  m.def("read_iges", &ReadIges, py::arg("fname"));
}

} // namespace splinepy::py
