#pragma once
// pybind11
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

// splinepy
#include "splinepy/fitting/fitting.hpp"

namespace splinepy::py {

namespace py = pybind11;

/// standalone version of fit_curve
void FitCurve(const py::array_t<double>& points,
              const int& degree,
              const int& num_control_points,
              const bool centripetal,
              const py::list& knot_vectors,
              py::dict& return_spline /* return */,
              double& residual /* return */);

/// @brief Interpolates curve through query points
/// @param points Query points
/// @param degree
/// @param centripetal
/// @param knot_vectors
/// @return py::dict
py::dict InterpolateCurve(py::array_t<double> points,
                          int degree,
                          bool centripetal,
                          py::list knot_vectors);

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
                          py::list knot_vectors);

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
                            bool centripetal);

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
                            bool centripetal);

} // namespace splinepy::py
