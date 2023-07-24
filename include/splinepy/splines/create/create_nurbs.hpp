#pragma once

#include <splinepy/splines/nurbs.hpp>
#include <splinepy/splines/splinepy_base.hpp>

namespace splinepy::splines::create {

/// dynamic creation of templated nurbs
std::shared_ptr<splinepy::splines::SplinepyBase>
CreateNurbs1(const int dim,
             const int* degrees,
             const std::vector<std::vector<double>>* knot_vectors,
             const double* control_points,
             const double* weights);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateNurbs2(const int dim,
             const int* degrees,
             const std::vector<std::vector<double>>* knot_vectors,
             const double* control_points,
             const double* weights);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateNurbs3(const int dim,
             const int* degrees,
             const std::vector<std::vector<double>>* knot_vectors,
             const double* control_points,
             const double* weights);

#ifdef SPLINEPY_MORE
std::shared_ptr<splinepy::splines::SplinepyBase>
CreateNurbs4(const int dim,
             const int* degrees,
             const std::vector<std::vector<double>>* knot_vectors,
             const double* control_points,
             const double* weights);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateNurbs5(const int dim,
             const int* degrees,
             const std::vector<std::vector<double>>* knot_vectors,
             const double* control_points,
             const double* weights);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateNurbs6(const int dim,
             const int* degrees,
             const std::vector<std::vector<double>>* knot_vectors,
             const double* control_points,
             const double* weights);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateNurbs7(const int dim,
             const int* degrees,
             const std::vector<std::vector<double>>* knot_vectors,
             const double* control_points,
             const double* weights);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateNurbs8(const int dim,
             const int* degrees,
             const std::vector<std::vector<double>>* knot_vectors,
             const double* control_points,
             const double* weights);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateNurbs9(const int dim,
             const int* degrees,
             const std::vector<std::vector<double>>* knot_vectors,
             const double* control_points,
             const double* weights);

std::shared_ptr<splinepy::splines::SplinepyBase>
CreateNurbs10(const int dim,
              const int* degrees,
              const std::vector<std::vector<double>>* knot_vectors,
              const double* control_points,
              const double* weights);
#endif

} // namespace splinepy::splines::create
