#include <sstream>

#include <BSplineLib/ParameterSpaces/knot_vector.hpp>

#include "splinepy/py/py_knot_insertion_matrix.hpp"
#include "splinepy/splines/helpers/extract.hpp"
#include "splinepy/utils/print.hpp"

/// @brief
namespace splinepy::py {

std::tuple<py::array_t<double>, std::vector<int>>
ComputeKnotInsertionMatrixAndKnotSpan(const py::array_t<double>& old_kv,
                                      const py::array_t<double>& new_kv,
                                      const int degree,
                                      const double& tolerance) {
  // Auxiliary access values
  const double* old_kv_ptr = static_cast<double*>(old_kv.request().ptr);
  const double* new_kv_ptr = static_cast<double*>(new_kv.request().ptr);
  const int n_cps_old{static_cast<int>(old_kv.size()) - degree - 1};
  const int n_cps_new{static_cast<int>(new_kv.size()) - degree - 1};
  // Create and initialize
  py::array_t<double> return_matrix(n_cps_new * n_cps_old);
  double* return_matrix_ptr = static_cast<double*>(return_matrix.request().ptr);
  for (int i{}; i < n_cps_new * n_cps_old; i++) {
    return_matrix_ptr[i] = 0.;
  }

  // Helper function, because I can't be bothered
  auto matrix = [&](const int& i, const int& j) -> double& {
    return return_matrix_ptr[i * n_cps_old + j];
  };

  // Helper function to create specific entries of the B-Spline Matrix on
  // the fly (as they are only required a single time)
  auto R_matrix = [&](const double& t, // Evaluation point
                      const int& mu,   // Offset with regards to local support
                      const int& deg,  // Matrix degree
                      const int& i,    // row
                      const int& j     // col
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
  std::vector<int> following_n_equal_knots(new_kv.size());
  std::vector<int> index_of_knot_span_base(new_kv.size());

  // As we assume closed knot-vectors (we can ignore the last entry)
  int current_old_id{static_cast<int>(old_kv.size() - 2)};
  int current_knot_span{static_cast<int>(old_kv.size() - 1)};

  // Prepare some values that will be required later (Ignore last value loop
  // begins at new_kv.size() - 2)
  for (int i{static_cast<int>(new_kv.size() - 1)}; i-- > 0;) {
    // Check knots against each other with a tolerance
    if (std::abs(old_kv_ptr[current_old_id] - new_kv_ptr[i]) < tolerance) {
      // Check if we entered another non-zero knot-span
      if (old_kv_ptr[current_old_id] < old_kv_ptr[current_old_id + 1]) {
        current_knot_span = current_old_id;
      }
      current_old_id--;
    } else {
      if (old_kv_ptr[current_old_id] > new_kv_ptr[i]) {
        // Beautify the output
        std::ostringstream old_s, new_s;
        new_s << "[" << new_kv_ptr[0];
        for (int i_kv{1}; i_kv < new_kv.size(); i_kv++) {
          new_s << ", " << new_kv_ptr[i];
        }
        new_s << "]";
        old_s << "[" << old_kv_ptr[0];
        for (int i_kv{1}; i_kv < old_kv.size(); i_kv++) {
          old_s << ", " << old_kv_ptr[i];
        }
        old_s << "]";
        splinepy::utils::PrintAndThrowError(
            "New knot spans is not subset of old knot-span. \nOld knot vector "
            ": ",
            old_s.str(),
            "\nNew knot vector : ",
            new_s.str());
      }
    }
    if (new_kv_ptr[i] < old_kv_ptr[current_knot_span]) {
      index_of_knot_span_base[i] = current_old_id;
    } else {
      index_of_knot_span_base[i] = current_knot_span;
    }
  }

  // Start assigning values to the matrix
  int offset{0};

  for (int i{0}; i < n_cps_new; i++) {
    // Get mu from previous knot_span
    const int mu{index_of_knot_span_base[i]};

    // Get to the first non-zero entry within the row
    int j{mu - degree};

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

    // Loop over all remaining matrices, use their sparsity to avoid unnecessary
    // computations and increase degree to be more efficient (no double matrix
    // computations)
    for (int d{2}; d <= degree; d++) {
      //  Loop over different matrices, last and first column of the d matrix
      //  must be treated differently, because they only have one entry, we loop
      //  backwards so we can avoid creating an additional auxiliary vectors
      matrix(i, j + d) =
          (R_matrix(new_kv_ptr[i + d], mu, d, d - 1, d) * matrix(i, j + d - 1));
      // For matrix multiplication we can go backwards to avoid creating an
      // additional aux vector and store solution of multiplication directly
      // into the solution matrix. A little harder to read, but makes the code
      // more efficient.
      for (int q{d - 1}; q > 0; q--) {
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
                           const int degree,
                           const double& tolerance) {
  return std::get<0>(
      ComputeKnotInsertionMatrixAndKnotSpan(old_kv, new_kv, degree, tolerance));
}

py::tuple ComputeGlobalKnotInsertionMatrix(
    const std::vector<py::array_t<double>>& old_kvs,
    const py::array_t<int>& degrees,
    const int parametric_dimension,
    const py::array_t<double>& new_knots,
    const double& tolerance) {

  // Create aliases
  const int* degrees_ptr = static_cast<int*>(degrees.request().ptr);
  const int n_para_dims = static_cast<int>(degrees.size());
  double* new_knots_ptr = static_cast<double*>(new_knots.request().ptr);

  // Make sure new knots are sorted
  std::sort(new_knots_ptr, new_knots_ptr + new_knots.size());

  // Precompute and initialize auxiliary values
  std::vector<int> n_cp_old(n_para_dims), n_cp_new(n_para_dims);
  int n_cps_old{1};
  for (int i{}; i < n_para_dims; i++) {
    // Knot vectors are stored as numpy arrays
    const py::array_t<double>& knot_vector_old = old_kvs[i];

    // Comput number of cps before and after insertion
    n_cp_old[i] = knot_vector_old.size() - degrees_ptr[i] - 1;
    n_cp_new[i] = n_cp_old[i];
    n_cps_old *= n_cp_old[i];
  }
  n_cp_new[parametric_dimension] += new_knots.size();
  const int n_cps_new = n_cps_old / n_cp_old[parametric_dimension]
                        * n_cp_new[parametric_dimension];

  // Compute dimensionwise (or line-) matrices. These matrices
  // are very small compared to higher dimensional spline matrices.
  const py::array_t<double>& knot_vector_old = old_kvs[parametric_dimension];
  double* knot_vector_old_ptr =
      static_cast<double*>(knot_vector_old.request().ptr);

  // Build new knot_vector
  py::array_t<double> knot_vector_new(knot_vector_old.size()
                                      + new_knots.size());
  double* knot_vector_new_ptr =
      static_cast<double*>(knot_vector_new.request().ptr);
  for (int i_old{}, j_new{}; i_old + j_new < knot_vector_new.size();) {
    if ((knot_vector_old_ptr[i_old] < new_knots_ptr[j_new])
        || (j_new >= new_knots.size())) {
      knot_vector_new_ptr[i_old + j_new] = knot_vector_old_ptr[i_old];
      i_old++;
    } else {
      knot_vector_new_ptr[i_old + j_new] = new_knots_ptr[j_new];
      j_new++;
    }
  }

  // Compute dimensionwise (or line-) matrices. These matrices
  // are very small compared to higher dimensional spline matrices.
  const auto results =
      ComputeKnotInsertionMatrixAndKnotSpan(knot_vector_old,
                                            knot_vector_new,
                                            degrees_ptr[parametric_dimension],
                                            tolerance);
  const py::array_t<double> dimensionwise_matrix{std::get<0>(results)};
  const double* dimensionwise_matrix_ptr =
      static_cast<double*>(dimensionwise_matrix.request().ptr);
  const std::vector<int> knot_span_ids{std::get<1>(results)};

  // Create auxiliary functions to access individual cps positions
  auto local_to_global_id = [&](const std::vector<int>& local_id,
                                const std::vector<int>& n_cps) -> int {
    int global_id{local_id[n_para_dims - 1]};
    for (int i{n_para_dims - 1}; i-- > 0;) {
      global_id *= n_cps[i + 0];
      global_id += local_id[i];
    }
    return global_id;
  };
  auto global_to_local_id =
      [&](const int& global_id,
          const std::vector<int>& n_cps) -> std::vector<int> {
    // Init
    std::vector<int> local_ids{};
    int copy_{global_id};
    local_ids.reserve(n_para_dims);
    for (int i{}; i < n_para_dims; i++) {
      const int r = copy_ % n_cps[i];
      local_ids.push_back(r);
      copy_ /= n_cps[i];
    }
    return local_ids;
  };

  // All matrices have specific size n_cps_new (post knot insertion) x
  // n_cps_old (prior to knot insertion), with a k-diagonal block structure,
  // where k is the degree along the parametric axis i. Some of the values
  // will be zero, if they do not contribute to one of the new inserted knots.

  // Initiate return values
  const int n_entries{(degrees_ptr[parametric_dimension] + 1) * n_cps_new};
  py::array_t<double> values(n_entries);
  py::array_t<int> rows(n_entries);
  py::array_t<int> cols(n_entries);

  // Retrieve pointers to data
  double* values_ptr = static_cast<double*>(values.request().ptr);
  int* rows_ptr = static_cast<int*>(rows.request().ptr);
  int* cols_ptr = static_cast<int*>(cols.request().ptr);

  // Assign values accordingly
  int counter{};
  for (int j_ctps{}; j_ctps < n_cps_new; ++j_ctps) {
    // j_ctps refers to the row associated to the current knot-span basis
    // n_cp_current
    std::vector<int> local_id_new = global_to_local_id(j_ctps, n_cp_new);
    const int i_row = local_id_new[parametric_dimension];
    // get knot-span
    const int knot_span_base_id =
        knot_span_ids[local_id_new[parametric_dimension]];

    for (int k_deg{knot_span_base_id - degrees_ptr[parametric_dimension]};
         k_deg <= knot_span_base_id;
         ++k_deg) {
      local_id_new[parametric_dimension] = k_deg;
      const int col_id = local_to_global_id(local_id_new, n_cp_old);
      // Assign values
      rows_ptr[counter] = static_cast<int>(j_ctps);
      cols_ptr[counter] = static_cast<int>(col_id);
      values_ptr[counter] =
          dimensionwise_matrix_ptr[i_row * n_cp_old[parametric_dimension]
                                   + k_deg];
      counter++;
    }
  }

  // Add contributions to list
  return py::make_tuple(
      py::make_tuple(values, py::make_tuple(rows, cols)),
      py::make_tuple(static_cast<int>(n_cps_new), static_cast<int>(n_cps_old)));
}

py::tuple BezierExtractionMatrices(const std::shared_ptr<PySpline>& spline,
                                   const double& tolerance) {
  // Create property copies
  const int n_para_dims = spline->Core()->SplinepyParaDim();
  std::vector<int> degrees(n_para_dims);
  std::vector<std::vector<double>> old_kvs;
  spline->Core()->SplinepyCurrentProperties(degrees.data(),
                                            &old_kvs,
                                            nullptr,
                                            nullptr);

  // Given a knot vector function provides a new knot vector in python format
  // that can represents the c^(-1) links
  auto create_new_knot_vector =
      [&](const py::array_t<double>& old_kv,
          const int& deg,
          py::array_t<double>& new_knots) -> py::array_t<double> {
    // allocate maximum required space for initialization
    const int old_kv_size = old_kv.size();
    const int max_kv_size = (old_kv_size - (deg + 1) * 2) * deg + (deg + 1) * 2;
    py::array_t<double> return_kv(max_kv_size);
    double* return_kv_ptr = static_cast<double*>(return_kv.request().ptr);
    double* old_kv_ptr = static_cast<double*>(old_kv.request().ptr);
    double* new_knots_ptr = static_cast<double*>(new_knots.request().ptr);
    int i_k_new{}, i_k_old{}, ref_count{};
    double c_knot{old_kv_ptr[0]};
    while (i_k_old < old_kv_size) {
      // If a jump to a new knot vector is a valid choice, we can update our
      // current knot reference
      if (ref_count > deg) {
        c_knot = old_kv_ptr[i_k_old];
        ref_count = 1;
      }

      if (std::abs(c_knot - old_kv_ptr[i_k_old]) < tolerance) {
        // Knot vector has not changed
        return_kv_ptr[i_k_new] = c_knot;
        ref_count++;
        i_k_old++;
        i_k_new++;
      } else {
        // Knot vector has changed
        return_kv_ptr[i_k_new] = c_knot;
        new_knots_ptr[i_k_new - i_k_old] = c_knot;
        i_k_new++;
        ref_count++;
      }
    }
    return_kv.resize({static_cast<int>(i_k_new)});
    new_knots.resize({static_cast<int>(i_k_new - i_k_old)});
    return return_kv;
  };

  // Initialize the return value
  py::list list_of_tuples{};
  // list of kv references - used to call ComputeGlobalKnotInsertionMatrix
  // we start with original kvs and this changes as we virtually keep inserting
  // knots each dimension
  std::vector<py::array_t<double>> tmp_kvs;
  tmp_kvs.reserve(n_para_dims);
  for (auto& kv : old_kvs) {
    // create views of old kvs as initial elements
    tmp_kvs.push_back(py::array_t<double>(kv.size(), kv.data(), py::none()));
  }

  // Here we compute the individual components that are required to initialize
  // the sparse matrices in the python frontend
  std::vector<int> n_patches_per_dimension{};
  n_patches_per_dimension.reserve(n_para_dims);
  // numpy view to call ComputeGlobalKnotInsertionMatrix
  py::array_t<int> np_degrees(degrees.size(), degrees.data(), py::none());
  for (int parametric_dimension{}; parametric_dimension < n_para_dims;
       parametric_dimension++) {
    const py::array_t<double>& tmp_kv = tmp_kvs[parametric_dimension];
    const int max_new_knots =
        (tmp_kv.size() - (degrees[parametric_dimension] + 1) * 2)
        * degrees[parametric_dimension];

    py::array_t<double> new_knots(max_new_knots);
    const auto& new_kv = create_new_knot_vector(tmp_kv,
                                                degrees[parametric_dimension],
                                                new_knots);
    n_patches_per_dimension.push_back(
        (new_kv.size() - degrees[parametric_dimension] - 2)
        / degrees[parametric_dimension]);

    // Add data to list
    list_of_tuples.append(ComputeGlobalKnotInsertionMatrix(tmp_kvs,
                                                           np_degrees,
                                                           parametric_dimension,
                                                           new_knots,
                                                           tolerance));
    tmp_kvs[parametric_dimension] = new_kv;
  }

  // to extract bezier ids, create knot multiplicity info.
  std::vector<std::vector<int>> knot_multiplicities;
  knot_multiplicities.reserve(n_para_dims);
  for (const py::array_t<double>& kv : tmp_kvs) {
    knot_multiplicities.emplace_back(
        bsplinelib::parameter_spaces::KnotVector::DetermineMultiplicities(
            static_cast<double*>(kv.request().ptr),
            kv.size(),
            tolerance));
  }

  // Lastly, compute the bezier extraction points, that correspond to the global
  // matrix.
  const auto& list_of_ids =
      splines::helpers::ExtractBezierPatchIDs(knot_multiplicities,
                                              degrees.data(),
                                              n_patches_per_dimension.data());
  const int n_patches = list_of_ids.size();
  const int n_ctps_per_patch = list_of_ids[0].size();
  py::array_t<int> bezier_ctps_ids(n_ctps_per_patch * n_patches);
  int* bezier_ctps_ids_ptr = static_cast<int*>(bezier_ctps_ids.request().ptr);
  for (int i_patch{}; i_patch < n_patches; i_patch++) {
    for (int j_id{}; j_id < n_ctps_per_patch; j_id++) {
      bezier_ctps_ids_ptr[i_patch * n_ctps_per_patch + j_id] =
          list_of_ids[i_patch][j_id];
    }
  }

  bezier_ctps_ids.resize(
      {static_cast<int>(n_patches), static_cast<int>(n_ctps_per_patch)});
  return py::make_tuple(bezier_ctps_ids, list_of_tuples);
}

// Provide function to add to module
void init_knot_insertion_matrix(py::module& m) {
  m.def("knot_insertion_matrix",
        &splinepy::py::ComputeKnotInsertionMatrix,
        py::arg("old_knot_vector"),
        py::arg("new_knot_vector"),
        py::arg("degree"),
        py::arg("tolerance"));
  m.def("global_knot_insertion_matrix",
        &splinepy::py::ComputeGlobalKnotInsertionMatrix,
        py::arg("old_knot_vectors"),
        py::arg("degrees"),
        py::arg("para_dim"),
        py::arg("new_knots"),
        py::arg("tolerance"));
  m.def("bezier_extraction_matrix",
        &splinepy::py::BezierExtractionMatrices,
        py::arg("spline"),
        py::arg("tolerance"));
}
} // namespace splinepy::py
