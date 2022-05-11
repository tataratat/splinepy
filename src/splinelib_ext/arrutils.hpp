#pragma once

#include <array>
#include <cmath>
#include <algorithm>

#include <Sources/Utilities/named_type.hpp>
#include <Sources/VectorSpaces/vector_space.hpp>
#include <Sources/ParameterSpaces/parameter_space.hpp>
  
using namespace splinelib::sources;  

/* elementwise subtraction */
template<typename T, int dim>
inline void ew_minus(const std::array<T, dim>& arr1,
                     const std::array<T, dim>& arr2,
                     std::array<T, dim>& out) {
  for (int i{0}; i < dim; i++) {
    out[i] = arr1[i] - arr2[i];
  }
}


/* elementwise subtraction overload for SplineLib Coordinate*/
template<typename T, int dim>
inline void ew_minus(
    const vector_spaces::VectorSpace<dim>::Coordinate_ arr1,
    const std::array<T, dim>& arr2,
    std::array<T, dim>& out) {
  for (int i{0}; i < dim; i++) {
    out[i] = arr1[splinelib::Index{i}].Get() - arr2[i];
  }
}


/* elementwise subtraction overload for SplineLib Coordinate*/
template<typename T, int dim>
inline void ew_minus(
    const std::array<T, dim>& arr1,
    const vector_spaces::VectorSpace<dim>::Coordinate_ arr2,
    std::array<T, dim>& out) {
  for (int i{0}; i < dim; i++) {
    out[i] = arr1[i] - arr2[splinelib::Index{i}].Get();
  }
}


/* elementwise subtraction overload for ptr */
template<typename T, int dim>
inline void ew_minus(
    const double* arr1,
    const vector_spaces::VectorSpace<dim>::Coordinate_ arr2,
    std::array<T, dim>& out) {

  for (int i{0}; i < dim; i++) {
    out[i] = arr1[i] - arr2[splinelib::Index{i}].Get();
  }

}


/* elementwise mean */
template<typename T, int dim>
inline void ew_mean(const std::array<T, dim>& arr1,
                    const std::array<T, dim>& arr2,
                    const std::array<T, dim>& out) {
  for (int i{0}; i < dim; i++) {
    out[i] = (arr1[i] + arr2[i]) * .5;
  }

}


template<typename T, int dim>
inline void norm2(std::array<T, dim>& in,
                  double& out) {
  double sqsum = 0.;
  for (int i{0}; i < dim; i++) {
    sqsum += in[i] * in[i];
  }

  out = std::sqrt(sqsum);
}

/* inplace operation for para coord clipping 
 *
 * Parameters
 * -----------
 * bounds: in
 * para_coord: out
 * clipped: out
 *   -1 -> minimum clip
 *    0 -> no clip
 *    1 -> maximum clip
 */
template<typename T, int para_dim>
inline void clip(
    std::array<std::array<double, para_dim>, 2>& bounds,
    parameter_spaces::ParameterSpace<para_dim>::ParameterCoordinate_& arr1,
    std::array<int, para_dim> clipped) {

  using PC =
    typename parameter_spaces::ParameterSpaces<para_dim>::ParameterCoordinate_;
  using SPC = PC::value_type;

  for (int i{0}; i < para_dim; i++) {
    // check max
    if (arr1[i] > bounds[1][i]) {
      clipped[i] = 1;
      arr1[i] = SPC{bounds[1][i]};
    } else if (arr1[i] < bounds[0][i]) {
      clipped[i] = -1;
      arr1[i] = SPC{bounds[0][i]};
    } else {
      clipped[i] = 0;
    }
  }
}


/* Gauss elimination with partial pivoting to find x from (A x = b)
 *
 * Parameters
 * -----------
 * A: in
 *   please excuse us, it will be modified
 * b: in
 *   please excuse us, it will be modified
 * x: out
 */
template<typename T, int para_dim>
inline void gauss_with_pivot(
    std::array<std::array<double, para_dim>, para_dim>& A,
    std::array<double, para_dim>& b,
    std::array<int, para_dim>& skipmask,
    std::array<double, para_dim>& x) {

  int maxrow;
  double maxval;
  
  // partial pivoting and forward reduction
  for (int i{0}; i < para_dim; i++) {
    // ignore clipped entries
    if (skipmask[i] != 0) continue;

    /* swap */
    // go through the rows and compare magnitude and mark it if needed
    maxrow = i;
    maxval = std::abs(A[i][i]);
    for (int j{i+1}; j < para_dim; j++) {
      curmax = std::abs(A[j][i]);
      if (curmax > maxval) {
        maxval = curmax;
        maxrow = j;
      }
    }
    // swap if needed. max entry row goes to i-th row.
    if (maxrow != i) {
      // swap matrix entries
      for (int j{0}; j < para_dim; j++) {
        std::swap(A[maxrow][j], A[i][j]);
      }
      // swap row-vec entries
      std::swap(b[maxrow], b[i]);
    }
    /* END swap */

    /* forward reduction */
    const double Aii = A[i][i];
    for (int j{i+1}; j < para_dim; j++) {
      const double redfac = A[j][i] / Aii;
      // matrix reduction
      for (int k{i+1}; k < para_dim; k++) {
        A[j][k] -= A[i][k] * redfac;
      }
      A[j][i] = 0.;

      // row-vec reduction0
      b[j] -= b[i] * redfac;
    }
    /* END forward reduction */

    // sneak in x initialization here.
    x[i] = 0.;
  }

  // back substitution
  double sum;
  for (int i{para_dim-1}; i >= 0; i--) {
    // ignore clipped entries
    if (skipmask[i] != 0) continue;

    sum = 0.;
    for (int j{i+1}; j < para_dim; j++) {
      sum += A[i][j] * x[j];
    }
    x[i] = (x[i] - sum) / A[i][i];
  }
}
