#pragma once

#include "lu_solve.hpp"
#include "utils.hpp"

namespace splinepy::fitting {
/**
 * @brief Interpolates a spline curve to points.
 *
 * @param points points array
 * @param num_points number of points
 * @param dim physical dimension
 * @param degree spline degree
 * @param num_control_points number of control points
 * @param centripetal centripetal parametrization
 * @param knot_vector knot vector
 * @param control_points control points array
 * @return double
 */
double FitCurve(const double* points,
                const int& num_points,
                const int& dim,
                const int& degree,
                const int& num_control_points,
                const bool centripetal,
                std::vector<double>& knot_vector,
                std::vector<double>& control_points);
/**
 * @brief Interpolates a spline surface to points.
 *
 * @param points rectangluar points array
 * @param num_points number of points
 * @param dim physical dimension
 * @param degree_u spline degree in u
 * @param degree_v spline degree in v
 * @param size_u number of control points in u
 * @param size_v number of control points in v
 * @param centripetal centripetal parametrization
 * @param knot_vector_u knot vector along u
 * @param knot_vector_v knot vector along v
 * @param control_points control points array
 */
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
/**
 * @brief Approximates a spline curve to points.
 *
 * @param points points
 * @param num_points number of points
 * @param dim physical dimension
 * @param degree spline degree
 * @param num_control_points number of control points
 * @param u_k parametrization
 * @param knot_vector knot vector
 * @param coefficient_matrix coefficient matrix of basis functions
 * @param control_points control points array
 * @return double
 */
double ApproximateCurve(const double* points,
                        const int& num_points,
                        const int& dim,
                        const int& degree,
                        const int& num_control_points,
                        std::vector<double>& u_k,
                        std::vector<double>& knot_vector,
                        std::vector<double>& coefficient_matrix,
                        std::vector<double>& control_points);

/**
 * @brief Approximates a spline surface to points shaped in a rectangular
 * grid along the parametric directions.
 *
 * @param points rectangluar points array
 * @param num_points_u number of points along u-direction
 * @param num_points_v number of points along v-direction
 * @param dim physical dimension
 * @param degree_u spline degree in u
 * @param degree_v spline degree in v
 * @param size_u number of control points in u
 * @param size_v number of control points in v
 * @param centripetal centripetal parametrization
 * @param knot_vector_u knot vector along u
 * @param knot_vector_v knot vector along v
 * @param control_points control points array
 */
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
