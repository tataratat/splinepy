#pragma once

#include <cmath>
#include <iterator>
#include <vector>

inline double square(double x) { return x * x; }

inline std::vector<double> ParametrizeCurve(const double* points,
                                            const int& num_points,
                                            const int& dim,
                                            const bool centripetal = true) {
  // Computes u_k
  int i, j;

  std::vector<double> chord_length;
  chord_length.assign(num_points + 1, 0.0);
  chord_length[num_points] = 1.0;

  double total_chord_length = 0.0;
  for (i = 1; i < num_points; i++) {
    for (j = 0; j < dim; j++) {
      chord_length[i] +=
          square(points[i * dim + j] - points[(i - 1) * dim + j]);
    }
    chord_length[i] = sqrt(chord_length[i]);
    if (centripetal) {
      chord_length[i] = sqrt(chord_length[i]);
    }
    total_chord_length += chord_length[i];
  }

  std::vector<double> u_k;
  u_k.assign(num_points, 0.0);

  for (i = 0; i < num_points; i++) {
    for (j = 0; j < i + 1; j++) {
      u_k[i] += chord_length[j];
    }
    u_k[i] /= total_chord_length;
  }

  return u_k;
}

inline void ParametrizeSurface(const double* points,
                               const int& num_points,
                               const int& dim,
                               const int& size_u,
                               const int& size_v,
                               const bool centripetal,
                               std::vector<double>& u_k,
                               std::vector<double>& v_l) {

  int i, u, v;
  double* pts_u = new double[size_u * dim];
  double* pts_v = new double[size_v * dim];
  std::vector<double> tmp_u_k{}, tmp_v_l{}, tmp_tmp_u_k, tmp_tmp_v_l;

  // Compute u_k
  u_k.assign(size_u, 0.0);

  // v - direction
  for (v = 0; v < size_v; v++) {
    for (u = 0; u < size_u; u++) {
      for (i = 0; i < dim; i++) {
        pts_u[u * dim + i] = points[(v + (size_v * u)) * dim + i];
      }
    }
    tmp_tmp_u_k = ParametrizeCurve(pts_u, size_u, dim, centripetal);
    std::move(tmp_tmp_u_k.begin(),
              tmp_tmp_u_k.end(),
              std::back_inserter(tmp_u_k));
  }

  // Average u - direction
  for (u = 0; u < size_u; u++) {
    for (v = 0; v < size_v; v++) {
      u_k[u] += tmp_u_k[u + (size_u * v)];
    }
    u_k[u] /= size_v;
  }

  // Compute u_k
  v_l.assign(size_v, 0.0);

  // u - direction
  for (u = 0; u < size_u; u++) {
    for (v = 0; v < size_v; v++) {
      for (i = 0; i < dim; i++) {
        pts_v[v * dim + i] = points[(v + (size_v * u)) * dim + i];
      }
    }
    tmp_tmp_v_l = ParametrizeCurve(pts_v, size_v, dim, centripetal);
    std::move(tmp_tmp_v_l.begin(),
              tmp_tmp_v_l.end(),
              std::back_inserter(tmp_v_l));
  }

  // Average v - direction
  for (v = 0; v < size_v; v++) {
    for (u = 0; u < size_u; u++) {
      v_l[v] += tmp_v_l[v + (size_v * u)];
    }
    v_l[v] /= size_u;
  }

  delete[] pts_u;
  delete[] pts_v;
}

inline std::vector<double> ComputeKnotVector(const int& degree,
                                             const int& num_points,
                                             const int& num_control_points,
                                             std::vector<double>& u_k) {

  int i, j;
  std::vector<double> knot_vector;
  knot_vector.assign(num_control_points + degree + 1, 0.0);

  // Distinguish between interpolation and approximation
  if (num_points == num_control_points) {

    // Equation 9.8 in the NURBS book
    for (i = degree + 1; i < num_points; i++) {
      for (j = i - degree; j < i; j++) {
        knot_vector[i] += u_k[j];
      }
      knot_vector[i] /= degree;
    }
  } else {
    // Equation 9,63 (note the index shift as num_control_points=n+1)
    const double d =
        ((double) num_points) / ((double) (num_control_points - degree));
    for (j = 1; j < num_control_points - degree; j++) {
      i = (int) (j * d);
      const double alpha = (j * d) - i;
      knot_vector[j + degree] = (1 - alpha) * u_k[i - 1] + alpha * u_k[i];
    }
  }

  for (i = num_control_points; i < num_control_points + degree + 1; i++) {
    knot_vector[i] = 1.0;
  }

  return knot_vector;
}

inline int FindSingleKnotSpan(const int& degree,
                              const std::vector<double>& knot_vector,
                              const int& num_control_points,
                              const double& knot) {

  int span = degree + 1;
  while (span < num_control_points && knot_vector[span] <= knot) {
    ++span;
  }

  return span - 1;
}

inline std::vector<double> BasisFunction(const int& degree,
                                         const std::vector<double>& knot_vector,
                                         const int& span,
                                         const double& knot) {

  std::vector<double> left, right, N;
  int i, j;
  left.assign(degree + 1, 0.0);
  right.assign(degree + 1, 0.0);
  N.assign(degree + 1, 1.0);

  double saved;
  double temp;
  for (i = 1; i < degree + 1; i++) {
    left[i] = knot - knot_vector[span + 1 - i];
    right[i] = knot_vector[span + i] - knot;

    saved = 0.0;
    for (j = 0; j < i; j++) {
      temp = N[j] / (right[j + 1] + left[i - j]);
      N[j] = saved + right[j + 1] * temp;
      saved = left[i - j] * temp;
    }
    N[i] = saved;
  }

  return N;
}

inline std::vector<double>
BuildCoefficientMatrix(const int& degree,
                       const std::vector<double>& knot_vector,
                       const std::vector<double>& u_k,
                       const int& num_points,
                       const int& num_control_points) {

  std::vector<double> coefficient_matrix;
  coefficient_matrix.assign(num_points * num_control_points, 0.0);
  std::vector<double> basis_function;
  int i, j, k, span;

  for (i = 0; i < num_points; i++) {
    span = FindSingleKnotSpan(degree, knot_vector, num_control_points, u_k[i]);
    basis_function = BasisFunction(degree, knot_vector, span, u_k[i]);
    k = 0;
    for (j = span - degree; j < span + 1; j++) {
      coefficient_matrix[i * num_control_points + j] = basis_function[k];
      k++;
    }
  }

  return coefficient_matrix;
}
