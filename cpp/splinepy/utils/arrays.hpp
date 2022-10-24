#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <type_traits>

namespace splinepy::utils {

/// Elementwise subtraction
template<typename T, typename InputArrayType, std::size_t dim>
inline void FirstMinusSecondEqualsThird(const InputArrayType& first,
                                        const T* second, /* c array */
                                        std::array<T, dim>& third) {
  // one exception for bezman's scalar splines
  if constexpr (std::is_scalar<InputArrayType>::value) {
    static_assert(dim == 1, "Minus with scalar is only for std::array<T, 1>");
    third[0] = first - second[0];
  } else {
    for (size_t i{}; i < dim; ++i) {
      // following line should raise error during compile time
      // if T != NameType::Type_
      third[i] = static_cast<T>(first[i]) - second[i];
    }
  }
}

/// Inplace version of Mean
/// Mainly to support bezman::Point types
template<typename ReturnArrayT, typename InputArrayT>
inline void
Mean(const InputArrayT& arr1, const InputArrayT& arr2, ReturnArrayT& out) {
  using ReturnValueType = typename ReturnArrayT::value_type;
  for (std::size_t i{}; i < out.size(); ++i) {
    out[i] = ReturnValueType{(arr1[i] + arr2[i]) * .5};
  }
}

/// Dot product - use SecondArrayType to support bezman's scalar splines
template<typename T1, typename SecondArrayT, std::size_t dim>
inline T1 Dot(const std::array<T1, dim>& arr1, const SecondArrayT& arr2) {
  T1 dotted{}; /* default value should be 0. */

  if constexpr (std::is_scalar<SecondArrayT>::value) {
    static_assert(dim == 1, "Dot with scalar is only for std::array<T, 1>");
    dotted = arr1[0] - arr2;
  } else {
    for (std::size_t i{0}; i < dim; ++i) {
      dotted += arr1[i] * static_cast<T1>(arr2[i]);
    }
  }

  return dotted;
}

/// AAt
template<typename T, std::size_t dim1, std::size_t dim2>
inline std::array<std::array<T, dim1>, dim1>
AAt(const std::array<std::array<T, dim2>, dim1>& arr1) {
  std::array<std::array<T, dim1>, dim1> out;

  for (std::size_t i{}; i < dim1; ++i) {
    for (std::size_t j{i}; j < dim1; ++j) {
      out[i][j] = Dot(arr1[i], arr1[j]);
      if (i != j)
        out[j][i] = out[i][j]; // fill symmetric part
    }
  }

  return out;
}

/// elementwise inplace addition
template<typename T1, typename T2, std::size_t dim>
inline void AddSecondToFirst(std::array<T1, dim>& arr1,
                             const std::array<T2, dim>& arr2) {
  for (std::size_t i{0}; i < dim; i++) {
    arr1[i] += T1{arr2[i]};
  }
}

/*!
 * Inplace operation for para coord clipping and saving clip info
 *
 * Parameters
 * -----------
 * @param[in] bounds
 * @param[out] para_coord
 * @param[out] clipped (-1) minimum clip; (0) no clip; (1) maximum clip
 */
template<typename T1, typename T2, std::size_t para_dim>
inline void Clip(const std::array<std::array<T1, para_dim>, 2>& bounds,
                 std::array<T2, para_dim>& para_coord,
                 std::array<int, para_dim>& clipped) {
  for (std::size_t i{0}; i < para_dim; i++) {
    // check max
    if (static_cast<T1>(para_coord[i]) > bounds[1][i]) {
      clipped[i] = 1;
      para_coord[i] = T2{bounds[1][i]};
      // check min
    } else if (static_cast<T1>(para_coord[i]) < bounds[0][i]) {
      clipped[i] = -1;
      para_coord[i] = T2{bounds[0][i]};
    } else {
      clipped[i] = 0;
    } // end if
  }   // end for
}

/// L2 Norm
template<typename T, std::size_t dim>
inline double NormL2(std::array<T, dim>& arr) {
  double returnval{};
  for (std::size_t i{}; i < dim; ++i) {
    returnval += arr[i] * arr[i];
  }
  return std::sqrt(returnval);
}

/* reorder by copying.
 * "who cares" approach
 * maybe faster. who knows
 *
 * Parameters
 * -----------
 * arr: inout  <- altered inplace
 * order: in
 */
template<typename T, typename IndexT, std::size_t dim>
void CopyReorder(std::array<T, dim>& arr, std::array<IndexT, dim>& order) {
  const auto copyarr = arr; // should copy.
  for (std::size_t i{0}; i < dim; i++) {
    arr[order[i]] = copyarr[i];
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
inline void
GaussWithPivot(std::array<std::array<double, para_dim>, para_dim>& A,
               std::array<double, para_dim>& b,
               std::array<int, para_dim>& skipmask,
               std::array<double, para_dim>& x) {
  std::size_t maxrow;
  double maxval, maybemax;
  std::array<int, para_dim> indexmap;
  std::iota(indexmap.begin(), indexmap.end(), 0);

  // partial pivoting and forward reduction
  for (std::size_t i{0}; i < para_dim; i++) {
    // sneak in x initialization here.
    // done first, so that skipping it will have no effect on update later
    x[i] = 0.;

    // ignore clipped entries
    if (skipmask[i] != 0)
      continue;

    /* swap */
    // go through the rows and compare magnitude and mark it if needed
    maxrow = i;
    maxval = std::abs(A[i][i]);
    for (std::size_t j{i + 1}; j < para_dim; j++) {
      maybemax = std::abs(A[j][i]);
      if (maybemax > maxval) {
        maxval = maybemax;
        maxrow = j;
      }
    }
    // swap if needed. max entry row goes to i-th row.
    if (maxrow != i) {
      // swap A's rows, entries of b, indexmap, and skipmask
      std::swap(A[maxrow], A[i]);
      std::swap(b[maxrow], b[i]);
      std::swap(indexmap[maxrow], indexmap[i]);
      std::swap(skipmask[maxrow], skipmask[i]);
    }
    /* END swap */

    /* forward reduction */
    const double Aii = A[i][i];
    for (std::size_t j{i + 1}; j < para_dim; j++) {
      const double redfac = A[j][i] / Aii;
      // matrix reduction
      for (std::size_t k{i + 1}; k < para_dim; k++) {
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
  // Looping backwards, stopping once overflow is reached, which corresponds to
  // (-1), i.e., such that the last value is 0
  for (std::size_t i{para_dim - 1}; i != static_cast<std::size_t>(-1); --i) {
    // ignore clipped entries
    if (skipmask[i] != 0)
      continue;

    sum = 0.;
    for (std::size_t j{i + 1}; j < para_dim; j++) {
      sum += A[i][j] * x[j];
    }
    x[i] = (b[i] - sum) / A[i][i];
  }

  // reorder
  CopyReorder(x, indexmap);
  CopyReorder(skipmask, indexmap);
}

template<typename T, std::size_t dim>
inline int NonZeros(std::array<T, dim>& arr) {
  int nonzeros{};
  for (const auto& a : arr) {
    if (a != 0)
      ++nonzeros;
  }
  return nonzeros;
}

} /* namespace splinepy::utils */
