#include <memory>
#include <vector>

#include "splinepy/py/py_spline_reader.hpp"
#include "splinepy/utils/print.hpp"

namespace splinepy::py {

/// @brief Convert BSpline to Python dictionary
/// @param bspline BSpline object
template<int para_dim, int dim>
py::dict
BSplineParser<para_dim, dim>::BSplineToDict(SplineEntry const& bspline) {
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
    ds_buf_ptr[i] = utilities::string_operations::ConvertToNumber<int>(degree);
    i++;
  }

  spline["degrees"] = p_degrees;

  // Unpack - knot vectors
  py::list p_knot_vectors;
  for (auto& knotvector : knot_vectors) {
    py::list p_kv;
    for (auto& knot : knotvector) {
      p_kv.append(utilities::string_operations::ConvertToNumber<double>(knot));
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

/// @brief BSpline parser class
/// @tparam para_dim Dimension of parametric space
/// @tparam dim Dimension of physical space
template<int para_dim, int dim>
py::dict NurbsParser<para_dim, dim>::NurbsToDict(SplineEntry const& nurbs) {
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
  OutputWeightedVectorSpace const& weighted_vector_space = std::get<1>(n_info);
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
    ds_buf_ptr[i] = utilities::string_operations::ConvertToNumber<int>(degree);
    i++;
  }

  spline["degrees"] = p_degrees;

  // Unpack - knot vectors
  py::list p_knot_vectors;
  for (auto& knotvector : knot_vectors) {
    py::list p_kv;
    for (auto& knot : knotvector) {
      p_kv.append(utilities::string_operations::ConvertToNumber<double>(knot));
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

py::list SplineReader::ReadIges(std::string fname) {
  int i, j;

  auto c_splines = input_output::iges::Read(fname);
  py::list dict_splines;

  for (auto& spline : c_splines) {
    if (spline->is_rational_) {
      dict_splines.append(ParseNurbs(spline));
    } else {
      dict_splines.append(ParseBSpline(spline));
    }
  }
  return dict_splines;
}

py::dict SplineReader::ParseBSpline(SplineEntry const& bspline) {

  // Get paradim and dim
  int const& para_dim = bspline->parametric_dimensionality_;
  int const& dim = bspline->dimensionality_;

  // Find appropriate BSplines
  switch (para_dim) {
  case 1:
    switch (dim) {
    case 1:
      return BSplineParser<1, 1>::BSplineToDict(bspline);
    case 2:
      return BSplineParser<1, 2>::BSplineToDict(bspline);
    case 3:
      return BSplineParser<1, 3>::BSplineToDict(bspline);
    default:
      splinepy::utils::PrintAndThrowError(
          "ParseBSpline() does not support physical dimension of ",
          dim,
          ". Only dimensions 1 to 3 are supported.");
      break;
    }
    break;
  case 2:
    switch (dim) {
    case 1:
      return BSplineParser<2, 1>::BSplineToDict(bspline);
    case 2:
      return BSplineParser<2, 2>::BSplineToDict(bspline);
    case 3:
      return BSplineParser<2, 3>::BSplineToDict(bspline);
    default:
      splinepy::utils::PrintAndThrowError(
          "ParseBSpline() does not support physical dimension of ",
          dim,
          ". Only dimensions 1 to 3 are supported.");
      break;
    }
    break;
  case 3:
    switch (dim) {
    case 1:
      return BSplineParser<3, 1>::BSplineToDict(bspline);
    case 2:
      return BSplineParser<3, 2>::BSplineToDict(bspline);
    case 3:
      return BSplineParser<3, 3>::BSplineToDict(bspline);
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
        ". Only dimensions 1 to 3 are supported.");
    break;
  }

  return py::dict{};
}

py::dict SplineReader::ParseNurbs(SplineEntry const& nurbs) {

  // Get paradim and dim
  int const& para_dim = nurbs->parametric_dimensionality_;
  int const& dim = nurbs->dimensionality_;

  // Find appropriate Nurbs
  switch (para_dim) {
  case 1:
    switch (dim) {
    case 1:
      return NurbsParser<1, 1>::NurbsToDict(nurbs);
    case 2:
      return NurbsParser<1, 2>::NurbsToDict(nurbs);
    case 3:
      return NurbsParser<1, 3>::NurbsToDict(nurbs);
    default:
      splinepy::utils::PrintAndThrowError(
          "ParseNurbs() does not support physical dimension of ",
          dim,
          ". Only dimensions 1 to 3 are supported.");
      break;
    }
    break;
  case 2:
    switch (dim) {
    case 1:
      return NurbsParser<2, 1>::NurbsToDict(nurbs);
    case 2:
      return NurbsParser<2, 2>::NurbsToDict(nurbs);
    case 3:
      return NurbsParser<2, 3>::NurbsToDict(nurbs);
    default:
      splinepy::utils::PrintAndThrowError(
          "ParseNurbs() does not support physical dimension of ",
          dim,
          ". Only dimensions 1 to 3 are supported.");
      break;
    }
    break;
  case 3:
    switch (dim) {
    case 1:
      return NurbsParser<3, 1>::NurbsToDict(nurbs);
    case 2:
      return NurbsParser<3, 2>::NurbsToDict(nurbs);
    case 3:
      return NurbsParser<3, 3>::NurbsToDict(nurbs);
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
        ". Only dimensions 1 to 3 are supported.");
    break;
  }

  return py::dict{};
}

///  @brief Adds spline reader. Keys are
/// ["knot_vectors", "control_points", "degrees"] (+ ["weights"])
/// @param m Python module
void init_spline_reader(py::module_& m) {
  m.def("read_iges", &SplineReader::ReadIges, py::arg("fname"));
}

} // namespace splinepy::py
