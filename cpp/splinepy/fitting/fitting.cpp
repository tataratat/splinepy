#include "fitting.hpp"

double FitCurve(double* points,
                int& num_points,
                int& dim,
                int& degree,
                int& num_control_points,
                bool centripetal,
                std::vector<double>& knot_vector,
                std::vector<double>& control_points) {

  std::vector<double> u_k, coefficient_matrix;

  u_k = ParametrizeCurve(points, num_points, dim, centripetal);

  // Check knot_vector dimensions and updated if required
  if (knot_vector.size()
      != static_cast<std::size_t>(num_control_points + degree + 1)) {
    knot_vector =
        ComputeKnotVector(degree, num_points, num_control_points, u_k);
  }

  coefficient_matrix = BuildCoefficientMatrix(degree,
                                              knot_vector,
                                              u_k,
                                              num_points,
                                              num_control_points);

  // Up until now there is no seperation between approximation and
  // interpolation.
  if (num_control_points == num_points) {
    control_points = LUSolve(coefficient_matrix, points, num_points, dim);
    return 0.0;
  } else {
    return ApproximateCurve(points,
                            num_points,
                            dim,
                            degree,
                            num_control_points,
                            u_k,
                            knot_vector,
                            coefficient_matrix,
                            control_points);
  }
}

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
                std::vector<double>& control_points) {

  std::vector<double> u_k, v_l, coefficient_matrix, tmp_result,
      tmp_control_points{};
  int u, v, i;
  double* pts_u = new double[size_u * dim];
  double* pts_v = new double[size_v * dim];

  // ParametrizeSurface
  ParametrizeSurface(points,
                     num_points,
                     dim,
                     size_u,
                     size_v,
                     centripetal,
                     u_k,
                     v_l);

  knot_vector_u = ComputeKnotVector(degree_u, size_u, size_u, u_k);

  knot_vector_v = ComputeKnotVector(degree_v, size_v, size_v, v_l);

  // u - direction global interpolation
  for (v = 0; v < size_v; v++) {
    for (u = 0; u < size_u; u++) {
      for (i = 0; i < dim; i++) {
        pts_u[u * dim + i] = points[(v + (size_v * u)) * dim + i];
      }
    }
    coefficient_matrix =
        BuildCoefficientMatrix(degree_u, knot_vector_u, u_k, size_u, size_u);
    tmp_result = LUSolve(coefficient_matrix, pts_u, size_u, dim);
    std::move(tmp_result.begin(),
              tmp_result.end(),
              std::back_inserter(tmp_control_points));
  }

  // v - direction global interpolation
  control_points.clear();
  for (u = 0; u < size_u; u++) {
    for (v = 0; v < size_v; v++) {
      for (i = 0; i < dim; i++) {
        pts_v[v * dim + i] = tmp_control_points[(u + (size_u * v)) * dim + i];
      }
    }
    coefficient_matrix =
        BuildCoefficientMatrix(degree_v, knot_vector_v, v_l, size_v, size_v);
    tmp_result = LUSolve(coefficient_matrix, pts_v, size_v, dim);
    std::move(tmp_result.begin(),
              tmp_result.end(),
              std::back_inserter(control_points));
  }

  delete[] pts_u;
  delete[] pts_v;
}

double ApproximateCurve(double* points,
                        int& num_points,
                        int& dim,
                        int& degree,
                        int& num_control_points,
                        std::vector<double>& u_k,
                        std::vector<double>& knot_vector,
                        std::vector<double>& coefficient_matrix,
                        std::vector<double>& control_points) {

  // Number of variable num_control_points (non-predetermined by P0 and Pm)
  int num_control_points_v = num_control_points - 2;

  // Calculate R_k
  double* R_points = new double[(num_points - 2) * dim];

  // equation 9.63
  for (int k{1}; k < num_points - 1; k++) {
    for (int i_dim{0}; i_dim < dim; i_dim++) {
      R_points[(k - 1) * dim + i_dim] =
          points[k * dim + i_dim]
          - coefficient_matrix[k * num_control_points] * points[i_dim]
          - coefficient_matrix[k * num_control_points + num_control_points - 1]
                * points[(num_points - 1) * dim + i_dim];
    }
  }

  // Assemble Least Square Problem
  // Reduced Matrix
  std::vector<double> LSmatrix;
  LSmatrix.assign(num_control_points_v * num_control_points_v, 0.0);

  // Reduced Vector NTR
  double* NtR = new double[num_control_points_v * dim];
  for (int i{}; i < num_control_points_v * dim; i++)
    NtR[i] = 0.0;

  // Calculate LS-System *symmetrize*
  for (int i{1}; i < num_control_points - 1; i++) {
    for (int k{1}; k < num_points - 1; k++) {
      // Assemble Reduced NTN
      for (int j{1}; j < num_control_points - 1; j++) {
        // Calculate NT N
        LSmatrix[(i - 1) * num_control_points_v + j - 1] +=
            coefficient_matrix[k * num_control_points + i]
            * coefficient_matrix[k * num_control_points + j];
      }

      // Calculate NtR
      for (int i_dim{0}; i_dim < dim; i_dim++) {
        NtR[(i - 1) * dim + i_dim] +=
            coefficient_matrix[k * num_control_points + i]
            * R_points[(k - 1) * dim + i_dim];
      }
    }
  }

  // Solve linear system with LU-decomposition
  control_points = LUSolve(LSmatrix, NtR, num_control_points_v, dim);

  // Add P0 at front and Pm at end
  control_points.resize(control_points.size() + dim * 2);
  for (int i{num_control_points_v * dim - 1}; i >= 0; i--) {
    control_points[i + dim] = control_points[i];
  }
  for (int i{0}; i < dim; i++) {
    control_points[i] = points[i];
    control_points[(num_control_points_v + 1) * dim + i] =
        points[(num_points - 1) * dim + i];
  }

  delete[] NtR;
  delete[] R_points;

  // Calculate Residual
  double residual = 0.0;
  for (int i{1}; i < num_points - 1; i++) {
    for (int j{0}; j < dim; j++) {
      double point_tmp{0.0};
      for (int k{0}; k < num_control_points; k++) {
        point_tmp += coefficient_matrix[i * num_control_points + k]
                     * control_points[k * dim + j];
      }

      residual +=
          (points[i * dim + j] - point_tmp) * (points[i * dim + j] - point_tmp);
    }
  }

  return residual;
}
