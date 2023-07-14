#pragma once

#include <splinepy/splines/rational_bezier.hpp>
#include <splinepy/splines/splinepy_base.hpp>

namespace splinepy::splines::create {

/// dynamic creation of templated rational
std::shared_ptr<splinepy::splines::SplinepyBase>
CreateRationalBezier(const int para_dim,
                     const int dim,
                     const int* degrees,
                     const double* control_points,
                     const double* weights);

} // namespace splinepy::splines::create
