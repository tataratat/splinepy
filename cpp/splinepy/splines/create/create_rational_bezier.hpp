#pragma once

#include <splinepy/splines/rational_bezier.hpp>
#include <splinepy/splines/splinepy_base.hpp>

namespace splinepy::splines::create {

/// dynamic creation of templated rational
std::shared_ptr<splinepy::splines::SplinepyBase>
CreateRationalBezier1(const int dim,
                      const int* degrees,
                      const double* control_points,
                      const double* weights);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateRationalBezier2(const int dim,
                      const int* degrees,
                      const double* control_points,
                      const double* weights);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateRationalBezier3(const int dim,
                      const int* degrees,
                      const double* control_points,
                      const double* weights);

#ifdef SPLINEPY_MORE
std::shared_ptr<splinepy::splines::SplinepyBase>
CreateRationalBezier4(const int dim,
                      const int* degrees,
                      const double* control_points,
                      const double* weights);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateRationalBezier5(const int dim,
                      const int* degrees,
                      const double* control_points,
                      const double* weights);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateRationalBezier6(const int dim,
                      const int* degrees,
                      const double* control_points,
                      const double* weights);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateRationalBezier7(const int dim,
                      const int* degrees,
                      const double* control_points,
                      const double* weights);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateRationalBezier8(const int dim,
                      const int* degrees,
                      const double* control_points,
                      const double* weights);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateRationalBezier9(const int dim,
                      const int* degrees,
                      const double* control_points,
                      const double* weights);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateRationalBezier10(const int dim,
                       const int* degrees,
                       const double* control_points,
                       const double* weights);
#endif

} // namespace splinepy::splines::create
