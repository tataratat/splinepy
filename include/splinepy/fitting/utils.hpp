#pragma once

#include <cmath>
#include <iterator>
#include <vector>

namespace splinepy::fitting {

inline double square(double const& x) { return x * x; }

std::vector<double> ParametrizeCurve(const double* points,
                                     const int& num_points,
                                     const int& dim,
                                     const bool centripetal = true);

void ParametrizeSurface(const double* points,
                        const int& num_points,
                        const int& dim,
                        const int& size_u,
                        const int& size_v,
                        const bool centripetal,
                        std::vector<double>& u_k,
                        std::vector<double>& v_l);

std::vector<double> ComputeKnotVector(const int& degree,
                                      const int& num_points,
                                      const int& num_control_points,
                                      std::vector<double>& u_k);

int FindSingleKnotSpan(const int& degree,
                       const std::vector<double>& knot_vector,
                       const int& num_control_points,
                       const double& knot);

std::vector<double> BasisFunction(const int& degree,
                                  const std::vector<double>& knot_vector,
                                  const int& span,
                                  const double& knot);

std::vector<double>
BuildCoefficientMatrix(const int& degree,
                       const std::vector<double>& knot_vector,
                       const std::vector<double>& u_k,
                       const int& num_points,
                       const int& num_control_points);

} // namespace splinepy::fitting
