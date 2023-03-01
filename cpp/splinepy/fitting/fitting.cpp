#include "fitting.hpp"

namespace splinepy::fitting {

double FitCurve(const double* points,
                const int& num_points,
                const int& dim,
                const int& degree,
                const int& num_control_points,
                const bool centripetal,
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
        pts_u[u * dim + i] = points[(u + (size_u * v)) * dim + i];
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
  control_points.assign(size_u * size_v * dim, 0.0);
  for (u = 0; u < size_u; u++) {
    for (v = 0; v < size_v; v++) {
      for (i = 0; i < dim; i++) {
        pts_v[v * dim + i] = tmp_control_points[(u + (size_u * v)) * dim + i];
      }
    }

    coefficient_matrix =
        BuildCoefficientMatrix(degree_v, knot_vector_v, v_l, size_v, size_v);
    tmp_result = LUSolve(coefficient_matrix, pts_v, size_v, dim);

    for (int v = 0; v < size_v; v++) {
      for (int i = 0; i < dim; i++) {
        control_points[(u + (size_u * v)) * dim + i] =
            tmp_result[(v * dim) + i];
      }
    }
  }

  delete[] pts_u;
  delete[] pts_v;
}

double ApproximateCurve(const double* points,
                        const int& num_points,
                        const int& dim,
                        const int& degree,
                        const int& num_control_points,
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
                        std::vector<double>& control_points) {

  std::vector<double> u_k, v_l, coefficient_matrix, tmp_result,
      tmp_control_points(size_u * num_points_v * dim),
      pts_u(num_points_u * dim), pts_v(num_points_v * dim);

  // Assign and compute the parameters u_k and v_l along both parametric
  // dimensions
  // Refer to NURBS Book Algorithm A9.3
  ParametrizeSurface(points,
                     num_points_u * num_points_v,
                     dim,
                     num_points_u,
                     num_points_v,
                     centripetal,
                     u_k,
                     v_l);

  knot_vector_u = ComputeKnotVector(degree_u, num_points_u, size_u, u_k);
  knot_vector_v = ComputeKnotVector(degree_v, num_points_v, size_v, v_l);

  // Build coefficient matrix containing the evaluated basis functions along
  // first parametric dimension
  // Refer to equation (9.66)
  coefficient_matrix = BuildCoefficientMatrix(degree_u,
                                              knot_vector_u,
                                              u_k,
                                              num_points_u,
                                              size_u);

  // Approximate temporary control point grid row-wise along the first
  // parametric direction as curves
  for (int v = 0; v < num_points_v; v++) {
    // Create pointer to avoid copying points
    const double* pointer_to_row = &points[v * num_points_u * dim];

    ApproximateCurve(pointer_to_row,
                     num_points_u,
                     dim,
                     degree_u,
                     size_u,
                     u_k,
                     knot_vector_u,
                     coefficient_matrix,
                     tmp_result);
    // Write result column-wise into the temporary control points
    for (int j = 0; j < size_u; j++) {
      for (int i = 0; i < dim; i++) {
        tmp_control_points[((num_points_v * j) + v) * dim + i] =
            tmp_result[(j * dim) + i];
      }
    }
  }

  coefficient_matrix.clear();
  coefficient_matrix = BuildCoefficientMatrix(degree_v,
                                              knot_vector_v,
                                              v_l,
                                              num_points_v,
                                              size_v);

  control_points.clear();
  control_points.assign(size_u * size_v * dim, 0.);

  // Interpolate each column of the temporary control grid along the second
  // parametric dimension
  for (int u = 0; u < size_u; u++) {

    // Create pointer to avoid copying points
    const double* pointer_to_column =
        &tmp_control_points[num_points_v * u * dim];

    ApproximateCurve(pointer_to_column,
                     num_points_v,
                     dim,
                     degree_v,
                     size_v,
                     v_l,
                     knot_vector_v,
                     coefficient_matrix,
                     tmp_result);

    // Write result column-wise into control points
    for (int v = 0; v < size_v; v++) {
      for (int i = 0; i < dim; i++) {
        control_points[(u + (size_u * v)) * dim + i] =
            tmp_result[(v * dim) + i];
      }
    }
  }
  tmp_result.clear();
}

} // namespace splinepy::fitting
