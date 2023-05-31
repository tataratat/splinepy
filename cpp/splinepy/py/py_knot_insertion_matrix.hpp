#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <splinepy/utils/print.hpp>

namespace splinepy::py {

namespace py = pybind11;

std::tuple<py::array_t<double>, std::vector<std::size_t>>
ComputeKnotInsertionMatrixAndKnotSpan(const py::array_t<double>& old_kv,
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
  for (std::size_t i{}; i < n_cps_new * n_cps_old; i++) {
    return_matrix_ptr[i] = 0.;
  }

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
  return std::make_tuple(return_matrix, index_of_knot_span_base);
}

py::array_t<double>
ComputeKnotInsertionMatrix(const py::array_t<double>& old_kv,
                           const py::array_t<double>& new_kv,
                           const std::size_t degree,
                           const double& tolerance) {
  return std::get<0>(
      ComputeKnotInsertionMatrixAndKnotSpan(old_kv, new_kv, degree, tolerance));
}

/**
 * @brief Compute the Global knot insertion matrix for all parametric dimensions
 * at once
 *
 * Currently we return a list of length (para_dim) containing a tuple of numpy
 * arrays, that can be used to instantiate scipy sparse matrices. This helps
 * avoid binding eigen as a seperate library.
 *
 * @todo: replace calls with pybind/eigen
 *
 * @param old_kvs list of arrays, representing the individual knot vectors of
 * the start spline
 * @param new_kvs list of arrays, representing the individual knot vectors of
 * the target spline
 * @param degrees degrees along all parametric dimensions
 * @param tolerance tolerance for idenifying individual knots
 * @return py::array_t<double> Conversion matrix
 */
py::list ComputeGlobalKnotInsertionMatrix(const py::list& old_kvs,
                                          const py::list& new_kvs,
                                          const py::array_t<int>& degrees,
                                          const double& tolerance) {

  // Create aliases
  const int* degrees_ptr = static_cast<int*>(degrees.request().ptr);
  const size_t n_para_dims = static_cast<size_t>(degrees.size());

  // Precompute and initialize auxiliary values. The dimensionwise matrices are
  // very small compared to higher dimensional spline matrices.
  std::vector<py::array_t<double>> dimensionwise_matrices{};
  std::vector<std::vector<std::size_t>> knot_span_ids{};
  std::vector<size_t> n_cp_old(n_para_dims), n_cp_new(n_para_dims);
  dimensionwise_matrices.reserve(n_para_dims);
  knot_span_ids.reserve(n_para_dims);
  std::size_t n_cps_old{1};
  for (std::size_t i{}; i < n_para_dims; i++) {
    // Knot vectors are stored as numpy arrays
    const py::array_t<double> knot_vector_old =
        py::cast<py::array_t<double>>(old_kvs[i]);
    const py::array_t<double> knot_vector_new =
        py::cast<py::array_t<double>>(new_kvs[i]);

    // Comput number of cps before and after insertion
    n_cp_old[i] = knot_vector_old.size() - degrees_ptr[i] - 1;
    n_cp_new[i] = knot_vector_new.size() - degrees_ptr[i] - 1;
    n_cps_old *= n_cp_old[i];

    // // Compute dimensionwise matrix and knot spans and assign to new values
    const auto results = ComputeKnotInsertionMatrixAndKnotSpan(knot_vector_old,
                                                               knot_vector_new,
                                                               degrees_ptr[i],
                                                               tolerance);
    dimensionwise_matrices.push_back(std::get<0>(results));
    knot_span_ids.push_back(std::get<1>(results));
  }

  // Create auxiliary functions to access individual cps positions
  auto local_to_global_id =
      [&](const std::vector<std::size_t>& local_id,
          const std::vector<std::size_t>& n_cps) -> std::size_t {
    std::size_t global_id{local_id[n_para_dims - 1]};
    for (std::size_t i{n_para_dims - 1}; i-- > 0;) {
      global_id *= n_cps[i + 0];
      global_id += local_id[i];
    }
    return global_id;
  };
  auto global_to_local_id =
      [&](const std::size_t& global_id,
          const std::vector<std::size_t>& n_cps) -> std::vector<std::size_t> {
    // Init
    std::vector<std::size_t> local_ids{};
    std::size_t copy_{global_id};
    local_ids.reserve(n_para_dims);
    for (std::size_t i{}; i < n_para_dims; i++) {
      const std::size_t r = copy_ % n_cps[i];
      local_ids.push_back(r);
      copy_ /= n_cps[i];
    }
    return local_ids;
  };

  // Loop over all parametric dimensions, each parametric dimension contributes
  // to a matrix.
  py::list return_value;
  // Keep track of current cps_dimensions
  std::vector<std::size_t> n_cp_current{n_cp_old};
  for (std::size_t i_para{}; i_para < n_para_dims; i_para++) {
    // All matrices have specific size n_cps_new (post knot insertion) x
    // n_cps_old (prior to knot insertion), with a k-diagonal block structure,
    // where k is the degree along the parametric axis i. Some of the values
    // will be zero, if they do not contribute to one of the new inserted knots.
    // We will update the new and old ctps size in every loop.
    std::size_t n_cps_new{n_cps_old / n_cp_old[i_para] * n_cp_new[i_para]};
    n_cp_current[i_para] = n_cp_new[i_para];

    // Auxiliary accessors
    const double* local_matrix =
        static_cast<double*>(dimensionwise_matrices[i_para].request().ptr);

    // Initiate return values
    const std::size_t n_entries{(degrees_ptr[i_para] + 1) * n_cps_new};
    py::array_t<double> values(n_entries);
    py::array_t<int> rows(n_entries);
    py::array_t<int> cols(n_entries);

    // Retrieve pointers to data
    double* values_ptr = static_cast<double*>(values.request().ptr);
    int* rows_ptr = static_cast<int*>(rows.request().ptr);
    int* cols_ptr = static_cast<int*>(cols.request().ptr);

    // Assign values accordingly
    std::size_t counter{};
    for (std::size_t j_ctps{}; j_ctps < n_cps_new; ++j_ctps) {
      // j_ctps refers to the row associated to the current knot-span basis
      // n_cp_current
      std::vector<std::size_t> local_id_new =
          global_to_local_id(j_ctps, n_cp_current);
      const std::size_t i_row = local_id_new[i_para];
      // get knot-span
      const std::size_t knot_span_base_id =
          knot_span_ids[i_para][local_id_new[i_para]];

      for (std::size_t k_deg{knot_span_base_id - degrees_ptr[i_para]};
           k_deg <= knot_span_base_id;
           ++k_deg) {
        local_id_new[i_para] = k_deg;
        const std::size_t col_id = local_to_global_id(local_id_new, n_cp_old);
        // Assign values
        rows_ptr[counter] = static_cast<int>(j_ctps);
        cols_ptr[counter] = static_cast<int>(col_id);
        values_ptr[counter] = local_matrix[i_row * n_cp_old[i_para] + k_deg];
        counter++;
      }
    }

    // Add contributions to list
    return_value.append(
        py::make_tuple(py::make_tuple(values, py::make_tuple(rows, cols)),
                       py::make_tuple(static_cast<int>(n_cps_new),
                                      static_cast<int>(n_cps_old))));

    // Update n_cp_old
    n_cps_old = n_cps_new;
    n_cp_old[i_para] = n_cp_new[i_para];
  }
  return return_value;
}
} // namespace splinepy::py
