#pragma once

#include "lu_solve.hpp"
#include "utils.hpp"

namespace splinepy::fitting {

double FitCurve(const double* points,
                const int& num_points,
                const int& dim,
                const int& degree,
                const int& num_control_points,
                const bool centripetal,
                std::vector<double>& knot_vector,
                std::vector<double>& control_points);

void FitSurface(const double* points,
                const int& num_points,
                const int& dim,
                const int& degree_u,
                const int& degree_v,
                const int& size_u,
                const int& size_v,
                const bool centripetal,
                std::vector<double>& knot_vector_u,
                std::vector<double>& knot_vector_v,
                std::vector<double>& control_points);

double ApproximateCurve(const double* points,
                        const int& num_points,
                        const int& dim,
                        const int& degree,
                        const int& num_control_points,
                        std::vector<double>& u_k,
                        std::vector<double>& knot_vector,
                        std::vector<double>& coefficient_matrix,
                        std::vector<double>& control_points);

void ApproximateSurface(const double* points,
                        const int& num_points_u,
                        const int& num_points_v,
                        const int& dim,
                        const int& degree_u,
                        const int& degree_v,
                        const int& size_u,
                        const int& size_v,
                        const bool centripetal,
                        std::vector<double>& knot_vector_u,
                        std::vector<double>& knot_vector_v,
                        std::vector<double>& control_points);
} // namespace splinepy::fitting
