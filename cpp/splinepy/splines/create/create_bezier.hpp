#pragma once

#include <splinepy/splines/bezier.hpp>
#include <splinepy/splines/rational_bezier.hpp>
#include <splinepy/splines/splinepy_base.hpp>

namespace splinepy::splines::create {

/// dynamic creation of templated Bezier
std::shared_ptr<splinepy::splines::SplinepyBase>
CreateBezier(const int para_dim,
             const int dim,
             const int* degrees,
             const double* control_points);

} // namespace splinepy::splines::create
