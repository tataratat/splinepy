#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <splinepy/utils/print.hpp>

namespace splinepy::py {

namespace py = pybind11;

py::array_t<double>
ComputeKnotInsertionMatrix(const py::array_t<double>& old_kv,
                           const py::array_t<double>& new_kv,
                           const std::size_t degree,
                           const double& tolerance) {

  // Auxiliary access values
  const double* old_kv_ptr = static_cast<double*>(old_kv.request().ptr);
  const double* new_kv_ptr = static_cast<double*>(new_kv.request().ptr);
  const std::size_t n_cps_old{old_kv.size() - degree - 1};
  const std::size_t n_cps_new{new_kv.size() - degree - 1};
  py::array_t<double> return_matrix(n_cps_new * n_cps_old);
  double* return_matrix_ptr = static_cast<double*>(return_matrix.request().ptr);

  // Helper function, because I can't be bothered
  auto matrix = [&](const std::size_t& i, const std::size_t& j) -> double& {
    return return_matrix_ptr[i * n_cps_old + j];
  };

  // Helper function to create specific entries of the B-Spline Matrix on
  // the fly (as they are only required a single time)
  auto R_matrix =
      [&](const double& t,        // Evaluation point
          const std::size_t& mu,  // Offset with regards to local support
          const std::size_t& deg, // Matrix degree
          const std::size_t& i,   // row
          const std::size_t& j    // col
          ) -> double {
    // Safety checks (Might be put into debugging block)
    if ((i >= deg) || j > deg) {
      splinepy::utils::PrintAndThrowError(
          "Requested index outside matrix definition.",
          "Requested query: \nt :",
          t,
          "\nmu : ",
          mu,
          "\ndeg : ",
          deg,
          "\ni : ",
          i,
          "\nj : ",
          j);
    }
    // Matrix is bi-diagonal i == j and i == j+1
    if ((i > j) || (j > i + 1)) {
      return 0.;
    }

    // Actual values
    if (j == i) {
      return (old_kv_ptr[mu + 1 + i] - t)
             / (old_kv_ptr[mu + 1 + i] - old_kv_ptr[mu + 1 + i - deg]);
    } else {
      return (t - old_kv_ptr[mu + 1 + i - deg])
             / (old_kv_ptr[mu + 1 + i] - old_kv_ptr[mu + 1 + i - deg]);
    }
  };

  // Another auxiliary values
  std::vector<std::size_t> following_n_equal_knots(new_kv.size());
  std::vector<std::size_t> index_of_knot_span_base(new_kv.size());

  // As we assume closed knot-vectors (we can ignore the last entry)
  std::size_t equal_index_count{0};
  std::size_t current_old_id{static_cast<std::size_t>(old_kv.size() - 2)};
  std::size_t current_knot_span{static_cast<std::size_t>(old_kv.size() - 1)};

  // Prepare some values that will be required later (Ignore last value loop
  // begins at new_kv.size() - 2)
  for (size_t i{static_cast<std::size_t>(new_kv.size() - 1)}; i-- > 0;) {
    // Check knots against each other with a tolerance
    if (std::abs(old_kv_ptr[current_old_id] - new_kv_ptr[i]) < tolerance) {
      equal_index_count++;
      // Check if we entered another non-zero knot-span
      if (old_kv_ptr[current_old_id] < old_kv_ptr[current_old_id + 1]) {
        current_knot_span = current_old_id;
      }
      current_old_id--;
    } else {
      if (old_kv_ptr[current_old_id] > new_kv_ptr[i]) {
        splinepy::utils::PrintAndThrowError(
            "New knot spans is not subset of old knot-span.");
      }
      equal_index_count = 0;
    }
    if (new_kv_ptr[i] < old_kv_ptr[current_knot_span]) {
      index_of_knot_span_base[i] = current_old_id;
    } else {
      index_of_knot_span_base[i] = current_knot_span;
    }
  }

  // Start assigning values to the matrix
  std::size_t offset{0};

  for (std::size_t i{0}; i < n_cps_new; i++) {
    // Get mu from previous knot_span
    const std::size_t mu{index_of_knot_span_base[i]};

    // Get to the first non-zero entry within the row
    std::size_t j{mu - degree};

    // Check if identity can be applied
    if (following_n_equal_knots[i + 1] > degree) {
      matrix(i, i - offset) = 1;
      // The degree == 0 case is also caught with this
      continue;
    }

    // Increase offset for every new knot
    if (following_n_equal_knots[i + 1] == 0) {
      offset++;
    }

    // As degree > 0 it must be at least 1 so we can pre-assign the first degree
    // values and write them into the matrix
    matrix(i, j) = R_matrix(new_kv_ptr[i + 1], mu, 1, 0, 0);
    matrix(i, j + 1) = R_matrix(new_kv_ptr[i + 1], mu, 1, 0, 1);

    // Loop over all remaining matrices, use their sparsity to avoid unecessary
    // computations and increase degree to be more efficient (no double matrix
    // computations)
    for (std::size_t d{2}; d <= degree; d++) {
      //  Loop over different matrices, last and first column of the d matrix
      //  must be treated differently, because they only have one entry, we loop
      //  backwards so we can avoid creating an additional auxiliary vectors
      matrix(i, j + d) =
          (R_matrix(new_kv_ptr[i + d], mu, d, d - 1, d) * matrix(i, j + d - 1));
      // For matrix multiplication we can go backwards to avoid creating an
      // additional aux vector and store solution of multiplication directly
      // into the solution matrix. A little harder to read, but makes the code
      // more efficient.
      for (std::size_t q{d - 1}; q > 0; q--) {
        matrix(i, j + q) =
            matrix(i, j + q - 1) * R_matrix(new_kv_ptr[i + d], mu, d, q - 1, q)
            + matrix(i, j + q) * R_matrix(new_kv_ptr[i + d], mu, d, q, q);
      }
      //  Assign first column last to complete matrix vector product for given
      //  line
      matrix(i, j) = R_matrix(new_kv_ptr[i + d], mu, d, 0, 0) * matrix(i, j);
    }
  }
  // assign shape and return result
  return_matrix.resize({(int) n_cps_new, (int) n_cps_old});
  return return_matrix;
}

} // namespace splinepy::py
