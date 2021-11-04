#include "utils.hpp"
#include "lu_solve.hpp"

double FitCurve(double* points,
                int& num_points,
                int& dim,
                int& degree,
                int& num_control_points,
                bool centripetal,
                std::vector<double>& knot_vector,
                std::vector<double>& control_points);

void FitSurface(double* points,
                int& num_points,
                int& dim,
                int& degree_u,
                int& degree_v,
                int& size_u,
                int& size_v,
                bool centripetal,
                std::vector<double>& knot_vector_u,
                std::vector<double>& knot_vector_v,
                std::vector<double>& control_points);

double ApproximateCurve(double* points,
                        int& num_points,
                        int& dim,
                        int& degree,
                        int& num_control_points,
                        std::vector<double>& u_k,
                        std::vector<double>& knot_vector,
                        std::vector<double>& coefficient_matrix,
                        std::vector<double>& control_points);
