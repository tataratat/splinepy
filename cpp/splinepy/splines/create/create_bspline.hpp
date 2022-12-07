#pragma once

#include <splinepy/splines/bspline.hpp>
#include <splinepy/splines/splinepy_base.hpp>

namespace splinepy::splines::create {

/// dynamic creation of templated BSpline
std::shared_ptr<splinepy::splines::SplinepyBase>
CreateBSpline1(const int dim,
               const int* degrees,
               const std::vector<std::vector<double>>* knot_vectors,
               const double* control_points);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateBSpline2(const int dim,
               const int* degrees,
               const std::vector<std::vector<double>>* knot_vectors,
               const double* control_points);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateBSpline3(const int dim,
               const int* degrees,
               const std::vector<std::vector<double>>* knot_vectors,
               const double* control_points);

#ifdef SPLINEPY_MORE
std::shared_ptr<splinepy::splines::SplinepyBase>
CreateBSpline4(const int dim,
               const int* degrees,
               const std::vector<std::vector<double>>* knot_vectors,
               const double* control_points);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateBSpline5(const int dim,
               const int* degrees,
               const std::vector<std::vector<double>>* knot_vectors,
               const double* control_points);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateBSpline6(const int dim,
               const int* degrees,
               const std::vector<std::vector<double>>* knot_vectors,
               const double* control_points);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateBSpline7(const int dim,
               const int* degrees,
               const std::vector<std::vector<double>>* knot_vectors,
               const double* control_points);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateBSpline8(const int dim,
               const int* degrees,
               const std::vector<std::vector<double>>* knot_vectors,
               const double* control_points);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateBSpline9(const int dim,
               const int* degrees,
               const std::vector<std::vector<double>>* knot_vectors,
               const double* control_points);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateBSpline10(const int dim,
                const int* degrees,
                const std::vector<std::vector<double>>* knot_vectors,
                const double* control_points);
#endif

} // namespace splinepy::splines::create
