#pragma once

#include <array>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iterator>

#include <Sources/Utilities/named_type.hpp>
#include <Sources/VectorSpaces/vector_space.hpp>
#include <Sources/ParameterSpaces/parameter_space.hpp>
  
using namespace splinelib::sources;  

/* elementwise subtraction */
template<typename T, std::size_t dim>
inline void ew_minus(const std::array<T, dim>& arr1,
                     const std::array<T, dim>& arr2,
                     std::array<T, dim>& out) {
  for (int i{0}; i < dim; i++) {
    out[i] = arr1[i] - arr2[i];
  }
}


/* elementwise subtraction overload for SplineLib Coordinate*/
template<typename T, std::size_t dim, typename VC /* VectorspaceCoordinate */>
inline void ew_minus(
    const VC& arr1,
    const std::array<T, dim>& arr2,
    std::array<T, dim>& out) {
  for (int i{0}; i < dim; i++) {
    out[i] = arr1[splinelib::Index{i}].Get() - arr2[i];
  }
}


/* elementwise subtraction overload for SplineLib Coordinate*/
template<typename T, std::size_t dim, typename VC>
inline void ew_minus(
    const std::array<T, dim>& arr1,
    const VC& arr2,
    std::array<T, dim>& out) {
  for (int i{0}; i < dim; i++) {
    out[i] = arr1[i] - arr2[splinelib::Index{i}].Get();
  }
}


/* elementwise subtraction overload for ptr */
template<typename T, std::size_t dim, typename VC>
inline void ew_minus(
    const double* arr1,
    const VC& arr2,
    std::array<T, dim>& out) {

  for (int i{0}; i < dim; i++) {
    out[i] = arr1[i] - arr2[splinelib::Index{i}].Get();
  }

}


/* elementwise mean */
template<typename T, std::size_t dim>
inline void ew_mean(const std::array<T, dim>& arr1,
                    const std::array<T, dim>& arr2,
                    std::array<T, dim>& out) {
  for (int i{0}; i < dim; i++) {
    out[i] = (arr1[i] + arr2[i]) * .5;
  }

}


/* elementwise addition */
template<typename T, std::size_t dim>
inline void ew_plus(const std::array<T, dim>& arr1,
                    const std::array<T, dim>& arr2,
                    std::array<T, dim>& out) {
  for (int i{0}; i < dim; i++) {
    out[i] = arr1[i] + arr2[i];
  }
}

/* inplace elementwise addition */
template<std::size_t para_dim, typename PC>
inline void ew_iplus(const std::array<double, para_dim>& arr1, PC& arr2) {

  using SPC = typename PC::value_type;

  int i{0};
  for (auto& a2 : arr2) {
    // TODO does one of them do less work?
    //a2 += SPC{arr1[i]};
    a2 = SPC{a2.Get() + arr1[i]};
  }
}


template<typename T, std::size_t dim>
inline void norm2(std::array<T, dim>& in,
                  double& out) {
  T sqsum{}; /* zero init */
  for (int i{0}; i < dim; i++) {
    sqsum += in[i] * in[i];
  }

  out = std::sqrt(sqsum);
}


template<typename T, std::size_t dim>
inline double norm2(std::array<T, dim>& in) {
  double returnval;
  norm2(in, returnval);
  return returnval;
}


/* inplace operation for para coord clipping 
 *
 * Parameters
 * -----------
 * bounds: in
 * para_coord: out
 * clipped: inout
 *   -1 -> minimum clip
 *    0 -> no clip
 *    1 -> maximum clip
 */
template<typename T, std::size_t para_dim, typename PC>
inline void clip(
    std::array<std::array<T, para_dim>, 2>& bounds,
    PC& arr1,
    std::array<int, para_dim> clipped) {

  using SPC = typename PC::value_type;

  for (int i{0}; i < para_dim; i++) {
    // check if it is already clipped
    if (clipped[i] != 0) continue;

    // check max
    if (arr1[i].Get() > bounds[1][i]) {
      clipped[i] = 1;
      arr1[i] = SPC{bounds[1][i]};
    } else if (arr1[i].Get() < bounds[0][i]) {
      clipped[i] = -1;
      arr1[i] = SPC{bounds[0][i]};
    } else { // I guess we won't reach here?
      clipped[i] = 0;
    }
  }
}


/* reorder. adapted from the world wide web.
 *
 * source:
 *  stackoverflow.com/questions/838384/reorder-vector-using-a-vector-of-indices
 * author: cgldr
 */
template<typename T, std::size_t dim>
void reorder(std::array<T, dim>& arr,
             std::array<int, dim>& order) {

  T ta;
  int i, j, k;
  /* reorder arr according to order */
  /* every move puts an element into place */
  /* time complexity is O(n) */
  for(i = 0; i < dim; i++){
    if(i != order[i]) {
      ta = arr[i];
      j = i;
      while(i != (k = order[j])) {
        arr[j] = arr[k];
        order[j] = j;
        j = k;
      }
      arr[j] = ta;
      order[j] = j;
    }
  }
}


/* Gauss elimination with partial pivoting to find x from (A x = b)
 *
 * Parameters
 * -----------
 * A: inout
 *   please excuse us, it will be modified
 * b: inout
 *   please excuse us, it will be modified
 * skipmask: in
 * x: out
 */
template<std::size_t para_dim>
inline void gauss_with_pivot(
    std::array<std::array<double, para_dim>, para_dim>& A,
    std::array<double, para_dim>& b,
    std::array<int, para_dim>& skipmask,
    std::array<double, para_dim>& x) {

  int maxrow;
  double maxval, maybemax;
  std::array<int, para_dim> indexmap;
  std::iota(indexmap.begin(), indexmap.end(), 0);
 
  // partial pivoting and forward reduction
  for (int i{0}; i < para_dim; i++) {
    // sneak in x initialization here.
    // done first, so that skipping it will have no effect on update later
    x[i] = 0.;

    // ignore clipped entries
    if (skipmask[i] != 0) continue;

    /* swap */
    // go through the rows and compare magnitude and mark it if needed
    maxrow = i;
    maxval = std::abs(A[i][i]);
    for (int j{i+1}; j < para_dim; j++) {
      maybemax = std::abs(A[j][i]);
      if (maybemax > maxval) {
        maxval = maybemax;
        maxrow = j;
      }
    }
    // swap if needed. max entry row goes to i-th row.
    if (maxrow != i) {
      // swap matrix entries
      // TODO: since it is nested, we can just do once?
      for (int j{0}; j < para_dim; j++) {
        std::swap(A[maxrow][j], A[i][j]);
      }
      // swap row-vec entries
      std::swap(b[maxrow], b[i]);
      // swap indexmap to return x correctly
      std::swap(indexmap[maxrow], indexmap[i]);
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

  // reorder
  reorder(x, indexmap);
}


template<typename T, std::size_t dim>
inline void dot(const std::array<T, dim>& arr1,
                const std::array<T, dim>& arr2,
                T& dotted) {
  dotted = T{}; /* default value should be 0. */
  for (int i{0}; i < dim; i++) {
    dotted += arr1[i] * arr2[i];
  }
}


/* maybe it is annoying sometimes to have dot as subroutine */
template <typename T, std::size_t dim>
inline T dot(const std::array<T, dim>& arr1,
             const std::array<T, dim>& arr2) {
  T returnval;
  dot(arr1, arr2, returnval);

  return returnval;
}


template<typename T, std::size_t dim, typename VC>
inline void dot(
    const std::array<T, dim>& arr1,
    const VC& arr2,
    T& dotted) {
  dotted = T{};
  for (int i{0}; i < dim; i++) {
    //dotted += arr1[i] * arr2[splinelib::Index{i}].Get();
    dotted += arr1[i] * arr2[i].Get();
  }
}


template<typename T, std::size_t dim, typename VC>
inline void dot(
    const VC& arr1,
    const std::array<T, dim>& arr2,
    T& dotted) {
  // switch place
  dot(arr2, arr1, dotted);
}


template<typename T, std::size_t dim, typename VC>
inline T dot(
    const std::array<T, dim>& arr1,
    const VC& arr2) {
  T returnval;
  dot(arr1, arr2, returnval);
  return returnval;
}


// makes `T dot` ambiguous. comment out 
//template<typename T, std::size_t dim, typename VC>
//inline T dot(
//    const VC& arr1,
//    const std::array<T, dim>& arr2) {
//  return dot(arr2, arr1);
//}


template<typename T, std::size_t dim1, std::size_t dim2>
inline void AAt(
    const std::array<std::array<T, dim2>, dim1>& arr1,
    std::array<std::array<T, dim1>, dim1>& out) {

  // use symmetry
  for (int i{0}; i < dim1; i++) {
    for (int j{0}; j < dim2; j++) {
      if (i > j) continue;
      out[i][j] = dot(arr1[i], arr1[j]);

      // fill symmetric part
      if (i != j) out[j][i] = out[i][j];
    }
  }

}


template<typename T, std::size_t dim>
inline int nonzeros(std::array<T, dim> arr) {
  int nnz{0};
  for (int i{0}; i < dim; i++) {
    if (arr[i] != 0) nnz++;
  }
  return nnz;
}
