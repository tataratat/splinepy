#pragma once
// pybind11
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

// splinepy
#include <splinepy/fitting/fitting.hpp>
#include <splinepy/utils/print.hpp>

namespace splinepy::py {

namespace py = pybind11;

/// standalone version of fit_curve
void FitCurve(const py::array_t<double>& points,
              const int& degree,
              const int& num_control_points,
              const bool centripetal,
              const py::list& knot_vectors,
              py::dict& return_spline /* return */,
              double& residual /* return */) {
  // Extract input array info.
  py::buffer_info p_buf = points.request();
  double* p_buf_ptr = static_cast<double*>(p_buf.ptr);

  // Prepare vars for interpolation
  int num_points = p_buf.shape[0];
  int curve_dim = p_buf.shape[1];

  // Update KnotVector (is optional)
  std::vector<double> knot_vector, control_points;
  for (auto& knotvector : knot_vectors) {
    for (auto& knot : knotvector) {
      // Empty knot_vector is indicator for updating
      knot_vector.push_back(knot.cast<double>());
    }
  }

  // write residual
  residual = splinepy::fitting::FitCurve(p_buf_ptr,
                                         num_points,
                                         curve_dim,
                                         degree,
                                         num_control_points,
                                         centripetal,
                                         knot_vector,
                                         control_points);

  // Write degree
  auto p_degrees = py::array_t<int>(1);
  py::buffer_info pd_buf = p_degrees.request();
  int* pd_buf_ptr = static_cast<int*>(pd_buf.ptr);
  pd_buf_ptr[0] = degree;

  // Write knot vector
  py::list p_knot_vectors;
  py::list kv;
  for (const auto& k : knot_vector) {
    kv.append(k);
  }
  p_knot_vectors.append(kv);

  // Write control points
  auto p_control_points = py::array_t<double>(control_points.size());
  py::buffer_info pc_buf = p_control_points.request();
  double* pc_buf_ptr = static_cast<double*>(pc_buf.ptr);

  for (int i_cps = 0; i_cps < num_control_points; i_cps++) {
    for (int j_dim = 0; j_dim < curve_dim; j_dim++) {
      pc_buf_ptr[i_cps * curve_dim + j_dim] =
          control_points[i_cps * curve_dim + j_dim];
    }
  }

  p_control_points.resize({num_control_points, curve_dim});

  // form return spline
  return_spline.attr("clear")();
  return_spline["degrees"] = p_degrees;
  return_spline["knot_vectors"] = p_knot_vectors;
  return_spline["control_points"] = p_control_points;
}

/// @brief Interpolates curve through query points
/// @param points Query points
/// @param degree
/// @param centripetal
/// @param knot_vectors
/// @return py::dict
py::dict InterpolateCurve(py::array_t<double> points,
                          int degree,
                          bool centripetal,
                          py::list knot_vectors) {
  py::buffer_info p_buf = points.request();
  int num_control_points = p_buf.shape[0];

  // prepare return values
  py::dict return_dict;
  double residual;

  // fit
  FitCurve(points,
           degree,
           num_control_points,
           centripetal,
           knot_vectors,
           return_dict,
           residual);

  return return_dict;
}

/// @brief Approximates curve based on query points
/// @param points Query points
/// @param degree
/// @param num_control_points
/// @param centripetal
/// @param knot_vectors
py::dict ApproximateCurve(py::array_t<double> points,
                          int degree,
                          int num_control_points,
                          bool centripetal,
                          py::list knot_vectors) {

  // prepare return values
  py::dict return_dict;
  double residual;
  FitCurve(points,
           degree,
           num_control_points,
           centripetal,
           knot_vectors,
           return_dict,
           residual);

  // append residual to it
  return_dict["residual"] = residual;

  return return_dict;
}

/// @brief Interpolates surface through query points
/// @param points Query points
/// @param size_u
/// @param size_v
/// @param degree_u
/// @param degree_v
/// @param centripetal
py::dict InterpolateSurface(py::array_t<double> points,
                            int size_u,
                            int size_v,
                            int degree_u,
                            int degree_v,
                            bool centripetal) {

  // Extract input array info.
  py::buffer_info p_buf = points.request();
  double* p_buf_ptr = static_cast<double*>(p_buf.ptr);

  // Prepare vars for interpolation
  int num_points = p_buf.shape[0];
  int surface_dim = p_buf.shape[1];
  if (surface_dim < 2) {
    splinepy::utils::PrintAndThrowError(
        "Query points should be at least 2D for surface fitting.");
  }
  std::vector<double> knot_vector_u, knot_vector_v, control_points;

  splinepy::fitting::FitSurface(p_buf_ptr,     //
                                num_points,    //
                                surface_dim,   //
                                degree_u,      //
                                degree_v,      //
                                size_u,        //
                                size_v,        //
                                centripetal,   //
                                knot_vector_u, //
                                knot_vector_v, //
                                control_points //
  );

  // Write degree
  auto p_degrees = py::array_t<int>(2);
  py::buffer_info pd_buf = p_degrees.request();
  int* pd_buf_ptr = static_cast<int*>(pd_buf.ptr);
  pd_buf_ptr[0] = degree_u;
  pd_buf_ptr[1] = degree_v;

  // Write knot vector
  py::list p_knot_vectors;
  py::list kv_u, kv_v;
  for (const auto& k : knot_vector_u) {
    kv_u.append(k);
  }
  for (const auto& k : knot_vector_v) {
    kv_v.append(k);
  }
  p_knot_vectors.append(kv_u);
  p_knot_vectors.append(kv_v);

  // Write control points
  auto p_control_points = py::array_t<double>(control_points.size());
  py::buffer_info pc_buf = p_control_points.request();
  double* pc_buf_ptr = static_cast<double*>(pc_buf.ptr);

  for (int i_cps = 0; i_cps < num_points; i_cps++) {
    for (int j_dim = 0; j_dim < surface_dim; j_dim++) {
      pc_buf_ptr[i_cps * surface_dim + j_dim] =
          control_points[i_cps * surface_dim + j_dim];
    }
  }

  p_control_points.resize({num_points, surface_dim});

  // prepare return
  py::dict dict_spline;
  dict_spline["degrees"] = p_degrees;
  dict_spline["knot_vectors"] = p_knot_vectors;
  dict_spline["control_points"] = p_control_points;

  return dict_spline;
}

/// @brief Approximates surface in the least-squares sense through query points
/// @param points The query points must form a rectangular grid along the x-
/// and y-axis
/// @param num_points_u The number of sampling points along the first
/// parametric direction. By default the first parametric direction is along
/// the cartesian x-axis, this can be adapted by reorganize.
/// @param num_points_v The number of sampling points along the second
/// parametric direction.
/// @param size_u Number of control points along first parametric direction
/// @param size_v Number of control points along second parametric direction
/// @param degree_u
/// @param degree_v
/// @param centripetal
py::dict ApproximateSurface(py::array_t<double> points,
                            int num_points_u,
                            int num_points_v,
                            int size_u,
                            int size_v,
                            int degree_u,
                            int degree_v,
                            bool centripetal) {

  // Extract input array info.
  py::buffer_info p_buf = points.request();
  double* p_buf_ptr = static_cast<double*>(p_buf.ptr);

  // Prepare vars for interpolation
  int surface_dim = p_buf.shape[1];
  if (surface_dim < 2) {
    splinepy::utils::PrintAndThrowError(
        "Query points should be at least 2D for surface fitting.");
  }
  std::vector<double> knot_vector_u, knot_vector_v, control_points;

  splinepy::fitting::ApproximateSurface(p_buf_ptr,     //
                                        num_points_u,  //
                                        num_points_v,  //
                                        surface_dim,   //
                                        degree_u,      //
                                        degree_v,      //
                                        size_u,        //
                                        size_v,        //
                                        centripetal,   //
                                        knot_vector_u, //
                                        knot_vector_v, //
                                        control_points //
  );
  // Write degree
  auto p_degrees = py::array_t<int>(2);
  py::buffer_info pd_buf = p_degrees.request();
  int* pd_buf_ptr = static_cast<int*>(pd_buf.ptr);
  pd_buf_ptr[0] = degree_u;
  pd_buf_ptr[1] = degree_v;

  // Write knot vector
  py::list p_knot_vectors;
  py::list kv_u, kv_v;
  for (const auto& k : knot_vector_u) {
    kv_u.append(k);
  }
  for (const auto& k : knot_vector_v) {
    kv_v.append(k);
  }
  p_knot_vectors.append(kv_u);
  p_knot_vectors.append(kv_v);

  // Write control points
  auto p_control_points = py::array_t<double>(control_points.size());
  py::buffer_info pc_buf = p_control_points.request();
  double* pc_buf_ptr = static_cast<double*>(pc_buf.ptr);

  for (int i_cps = 0; i_cps < (size_u * size_v); i_cps++) {
    for (int j_dim = 0; j_dim < surface_dim; j_dim++) {
      pc_buf_ptr[i_cps * surface_dim + j_dim] =
          control_points[i_cps * surface_dim + j_dim];
    }
  }

  p_control_points.resize({size_u * size_v, surface_dim});

  // prepare return
  py::dict dict_spline;
  dict_spline["degrees"] = p_degrees;
  dict_spline["knot_vectors"] = p_knot_vectors;
  dict_spline["control_points"] = p_control_points;

  return dict_spline;
}

/// @brief Functions that return fitted BSpline as dict
/// @param m
inline void add_fitting(py::module& m) {
  // Functions that return fitted bspline as dict.
  m.def("interpolate_curve",
        &splinepy::py::InterpolateCurve,
        py::arg("points"),
        py::arg("degree"),
        py::arg("centripetal"),
        py::arg("knot_vector"));
  m.def("approximate_curve",
        &splinepy::py::ApproximateCurve,
        py::arg("points"),
        py::arg("degree"),
        py::arg("n_control_points"),
        py::arg("centripetal"),
        py::arg("knot_vector"));
  m.def("interpolate_surface",
        &splinepy::py::InterpolateSurface,
        py::arg("points"),
        py::arg("size_u"),
        py::arg("size_v"),
        py::arg("degree_u"),
        py::arg("degree_v"),
        py::arg("centripetal"));
  m.def("approximate_surface",
        &splinepy::py::ApproximateSurface,
        py::arg("points"),
        py::arg("num_points_u"),
        py::arg("num_points_v"),
        py::arg("size_u"),
        py::arg("size_v"),
        py::arg("degree_u"),
        py::arg("degree_v"),
        py::arg("centripetal"));
}

} // namespace splinepy::py
