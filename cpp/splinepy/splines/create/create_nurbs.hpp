#pragma once

#include <splinepy/splines/nurbs.hpp>
#include <splinepy/splines/splinepy_base.hpp>

namespace splinepy::splines::create {

/// dynamic creation of templated nurbs
std::shared_ptr<splinepy::splines::SplinepyBase>
CreateNurbs(const int para_dim,
            const int dim,
            const int* degrees,
            const std::vector<std::vector<double>>* knot_vectors,
            const double* control_points,
            const double* weights);

} // namespace splinepy::splines::create
