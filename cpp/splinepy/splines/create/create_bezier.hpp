#pragma once

#include <splinepy/splines/bezier.hpp>
#include <splinepy/splines/rational_bezier.hpp>
#include <splinepy/splines/splinepy_base.hpp>

namespace splinepy::splines::create {

/// dynamic creation of templated Bezier
std::shared_ptr<splinepy::splines::SplinepyBase>
CreateBezier1(const int dim, const int* degrees, const double* control_points);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateBezier2(const int dim, const int* degrees, const double* control_points);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateBezier3(const int dim, const int* degrees, const double* control_points);

#ifdef SPLINEPY_MORE
std::shared_ptr<splinepy::splines::SplinepyBase>
CreateBezier4(const int dim, const int* degrees, const double* control_points);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateBezier5(const int dim, const int* degrees, const double* control_points);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateBezier6(const int dim, const int* degrees, const double* control_points);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateBezier7(const int dim, const int* degrees, const double* control_points);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateBezier8(const int dim, const int* degrees, const double* control_points);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateBezier9(const int dim, const int* degrees, const double* control_points);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateBezier10(const int dim, const int* degrees, const double* control_points);
#endif

} // namespace splinepy::splines::create
