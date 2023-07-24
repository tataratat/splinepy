#pragma once

#include <cmath>
#include <vector>

namespace splinepy::fitting {

void Doolittle(const std::vector<double>& matrix,
               std::vector<double>& l,
               std::vector<double>& u);

void ForwardSubstitution(const std::vector<double>& l,
                         const std::vector<double>& b,
                         std::vector<double>& y,
                         const int& num_b);

void BackwardSubstitution(const std::vector<double>& u,
                          const std::vector<double>& y,
                          std::vector<double>& x,
                          const int& num_b);

std::vector<double> LUSolve(const std::vector<double>& matrix,
                            const double* b,
                            const int& num_b,
                            const int& b_dim);

} // namespace splinepy::fitting
