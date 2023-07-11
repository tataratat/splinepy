#pragma once

#include <splinepy/splines/bspline.hpp>
#include <splinepy/splines/splinepy_base.hpp>

namespace splinepy::splines::create {

/// dynamic creation of templated BSpline
std::shared_ptr<splinepy::splines::SplinepyBase>
CreateBSpline(const int para_dim,
              const int dim,
              const int* degrees,
              const std::vector<std::vector<double>>* knot_vectors,
              const double* control_points);
} // namespace splinepy::splines::create
