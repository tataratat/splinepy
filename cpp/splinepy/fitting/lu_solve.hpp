#pragma once

#include <cmath>
#include <vector>

inline void Doolittle(std::vector<double>& matrix,
                      std::vector<double>& l,
                      std::vector<double>& u) {

  int i, j, k, num_rows = sqrt(matrix.size());
  for (i = 0; i < num_rows; i++) {
    for (k = 0; k < num_rows; k++) {
      for (j = 0; j < i; j++) {
        u[i * num_rows + k] -= l[i * num_rows + j] * u[j * num_rows + k];
      }
      u[i * num_rows + k] += matrix[i * num_rows + k];

      if (i == k) {
        l[i * num_rows + i] = 1.0;
      } else {
        for (j = 0; j < i; j++) {
          l[k * num_rows + i] -= l[k * num_rows + j] * u[j * num_rows + i];
        }
        l[k * num_rows + i] += matrix[k * num_rows + i];

        if (u[i * num_rows + i] == 0.0) {
          l[k * num_rows + i] = 0.0;
        } else {
          l[k * num_rows + i] /= u[i * num_rows + i];
        }
      }
    }
  }
}

inline void ForwardSubstitution(std::vector<double>& l,
                                std::vector<double>& b,
                                std::vector<double>& y,
                                int& num_b) {

  int i, j;
  y[0] = b[0] / l[0];
  for (i = 1; i < num_b; i++) {
    for (j = 0; j < i; j++) {
      y[i] -= l[i * num_b + j] * y[j];
    }
    y[i] += b[i];
    y[i] /= l[i * num_b + i];
  }
}

inline void BackwardSubstitution(std::vector<double>& u,
                                 std::vector<double>& y,
                                 std::vector<double>& x,
                                 int& num_b) {

  int i, j;
  x[num_b - 1] = y[num_b - 1]
                 / u[num_b * num_b - 1]; // equiv. `u`'s last elem ( u[-1,-1] )
  for (i = num_b - 2; i > -1; i--) {
    for (j = i; j < num_b; j++) {
      x[i] -= u[i * num_b + j] * x[j];
    }
    x[i] += y[i];
    x[i] /= u[i * num_b + i];
  }
}

inline std::vector<double>
LUSolve(std::vector<double>& matrix, double* b, int& num_b, int& b_dim) {

  std::vector<double> x, l, u, tmp_b, y, tmp_x;
  int i, j;

  x.assign(num_b * b_dim, 0.0);
  l.assign(matrix.size(), 0.0);
  u.assign(matrix.size(), 0.0);

  tmp_b.assign(num_b, 0.0);

  Doolittle(matrix, l, u);

  for (i = 0; i < b_dim; i++) {
    for (j = 0; j < num_b; j++) {
      tmp_b[j] = b[j * b_dim + i];
    }
    y.assign(num_b, 0.0);
    ForwardSubstitution(l, tmp_b, y, num_b);
    tmp_x.assign(num_b, 0.0);
    BackwardSubstitution(u, y, tmp_x, num_b);
    for (j = 0; j < num_b; j++) {
      x[j * b_dim + i] = tmp_x[j];
    }
  }

  return x;
}
