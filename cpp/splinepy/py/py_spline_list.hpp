#pragma once

#include <cstdlib>
#include <string>

// pybind
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

#include <splinepy/py/py_spline.hpp>
#include <splinepy/py/py_spline_extensions.hpp>
#include <splinepy/splines/helpers/scalar_type_wrapper.hpp>
#include <splinepy/utils/default_initialization_allocator.hpp>
#include <splinepy/utils/nthreads.hpp>

namespace splinepy::py {

namespace py = pybind11;

// alias
using PySplineList = std::vector<std::shared_ptr<PySpline>>;
using IntVector = splinepy::utils::DefaultInitializationVector<int>;
using DoubleVector = splinepy::utils::DefaultInitializationVector<double>;

/// @brief raises if elements' para_dim and dims aren't equal
/// @param splist
/// @param nthreads
inline void RaiseIfElementsHaveUnequalParaDimOrDim(const PySplineList& splist,
                                                   const int nthreads) {
  // use first spline as guide line
  const int& para_dim = splist[0]->para_dim_;
  const int& dim = splist[0]->dim_;

  std::vector<IntVector> mismatches(nthreads);

  auto check_dims_step = [&](int begin, int total_) {
    for (int i{begin}; i < total_; i += nthreads) {
      if (!CheckParaDimAndDim(*splist[i], para_dim, dim, false)) {
        mismatches[begin].push_back(i);
      }
    }
  };

  splinepy::utils::NThreadExecution(check_dims_step,
                                    static_cast<int>(splist.size()),
                                    nthreads,
                                    splinepy::utils::NThreadQueryType::Step);
  // prepare error prints
  IntVector all_mismatches{};
  for (const auto& m : mismatches) {
    const auto m_size = m.size();
    if (m_size != 0) {
      all_mismatches.reserve(all_mismatches.size() + m_size);
      all_mismatches.insert(all_mismatches.end(), m.begin(), m.end());
    }
  }

  if (all_mismatches.size() == 0) {
    return;
  }

  std::string mismatch_ids{};
  for (const auto& am : all_mismatches) {
    mismatch_ids += std::to_string(am);
    mismatch_ids += ", ";
  }
  splinepy::utils::PrintAndThrowError(
      "Splines in following entries has mismatching para_dim or dim compared "
      "to the first entry (",
      mismatch_ids,
      ")");
}

/// evaluates splines of same para_dim and dim.
inline py::array_t<double> ListEvaluate(const PySplineList& splist,
                                        const py::array_t<double>& queries,
                                        const int nthreads,
                                        const bool check_dims) {
  // use first spline as dimension guide line
  const int& para_dim = splist[0]->para_dim_;
  const int& dim = splist[0]->dim_;

  // query dim check
  CheckPyArrayShape(queries, {-1, para_dim}, true);

  // check dims if wanted. Else it will assume that every spline has same
  // dimension as the first entry's
  if (check_dims) {
    RaiseIfElementsHaveUnequalParaDimOrDim(splist, nthreads);
  }

  // prepare input and output
  double* queries_ptr = static_cast<double*>(queries.request().ptr);
  const int n_splines = splist.size();
  const int n_queries = queries.shape(0);
  const int n_total = n_splines * n_queries;
  py::array_t<double> evaluated({n_total, dim});
  double* evaluated_ptr = static_cast<double*>(evaluated.request().ptr);

  // each thread evaluates similar amount of queries from each spline
  auto evaluate_step = [&](int begin, int total_) {
    for (int i{begin}; i < total_; i += nthreads) {
      const auto [i_spline, i_query] = std::div(i, n_queries);
      splist[i_spline]->Core()->SplinepyEvaluate(
          &queries_ptr[i_query * para_dim],
          &evaluated_ptr[(i_spline * n_queries + i_query) * dim]);
    }
  };

  // exe
  splinepy::utils::NThreadExecution(evaluate_step,
                                    n_total,
                                    nthreads,
                                    splinepy::utils::NThreadQueryType::Step);

  return evaluated;
}

inline py::array_t<double> ListDerivative(const PySplineList& splist,
                                          const py::array_t<double>& queries,
                                          const py::array_t<int>& orders,
                                          const int nthreads) {
  return py::array_t<double>();
}

/// Samples equal resoltions for each para dim
inline py::array_t<double> ListSample(const PySplineList& splist,
                                      const int resolution,
                                      const int nthreads,
                                      const bool same_parametric_bounds,
                                      const bool check_dims) {
  // use first spline as dimension guide line
  const auto& first_spline = *splist[0];
  const int& para_dim = first_spline.para_dim_;
  const int& dim = first_spline.dim_;

  // dim check first
  if (check_dims) {
    RaiseIfElementsHaveUnequalParaDimOrDim(splist, nthreads);
  }

  // n_queries, and n_splines;
  int n_queries{1};
  const int n_splines = splist.size();

  // prepare resolutions
  IntVector resolutions_vector(para_dim);
  int* resolutions = resolutions_vector.data();
  for (int i{}; i < para_dim; ++i) {
    resolutions[i] = resolution;
    n_queries *= resolution;
  }

  // prepare input /  output
  const int n_total = n_splines * n_queries;
  py::array_t<double> sampled({n_total, dim});
  double* sampled_ptr = static_cast<double*>(sampled.request().ptr);

  // if you know all the queries have same parametric bounds
  // you don't need to re-compute queries
  if (same_parametric_bounds) {

    // get para bounds
    DoubleVector para_bounds_vector(2 * para_dim);
    double* para_bounds = para_bounds_vector.data();
    first_spline.Core()->SplinepyParametricBounds(para_bounds);

    // prepare queries
    DoubleVector queries_vector(n_queries * para_dim);
    double* queries = queries_vector.data();

    // use grid point generator to fill queries
    splinepy::utils::CStyleArrayPointerGridPoints gp_generator(para_dim,
                                                               para_bounds,
                                                               resolutions);
    gp_generator.Fill(queries);

    // create lambda for nthread exe
    auto sample_same_bounds_step = [&](int begin, int total_) {
      for (int i{begin}; i < total_; i += nthreads) {
        const auto [i_spline, i_query] = std::div(i, n_queries);
        splist[i_spline]->Core()->SplinepyEvaluate(
            &queries[i_query * para_dim],
            &sampled_ptr[(i_spline * n_queries + i_query) * dim]);
      }
    };

    splinepy::utils::NThreadExecution(sample_same_bounds_step,
                                      n_total,
                                      nthreads,
                                      splinepy::utils::NThreadQueryType::Step);

  } else {
    // here, we will execute 2 times:
    //   first, to create grid point helpers for each spline
    //   second, to sample

    // create a container to hold grid point helper.
    splinepy::utils::DefaultInitializationVector<
        splinepy::utils::CStyleArrayPointerGridPoints>
        grid_points(n_splines);

    // create grid_points
    auto create_grid_points = [&](int begin, int end) {
      DoubleVector para_bounds_vector(2 * para_dim);
      double* para_bounds = para_bounds_vector.data();

      for (int i{begin}; i < end; ++i) {
        // get para_bounds
        splist[i]->Core()->SplinepyParametricBounds(para_bounds);
        // setup grid points helper
        grid_points[i].SetUp(para_dim, para_bounds, resolutions);
      }
    };

    // pre compute entries -> this one is a chunk query
    splinepy::utils::NThreadExecution(create_grid_points, n_splines, nthreads);

    // similar to the one with same_parametric_bounds, except it computes query
    // on the fly
    auto sample_step = [&](int begin, int total_) {
      // each thread needs just one query array
      DoubleVector thread_query_vector(para_dim);
      double* thread_query = thread_query_vector.data();

      for (int i{begin}; i < total_; i += nthreads) {
        const auto [i_spline, i_query] = std::div(i, n_queries);
        const auto& gp_helper = grid_points[i_spline];
        gp_helper.IdToGridPoint(i_query, thread_query);
        splist[i_spline]->Core()->SplinepyEvaluate(
            thread_query,
            &sampled_ptr[(i_spline * n_queries + i_query) * dim]);
      }
    };

    // exe - this one is step
    splinepy::utils::NThreadExecution(sample_step,
                                      n_total,
                                      nthreads,
                                      splinepy::utils::NThreadQueryType::Step);
  }

  return sampled;
}

// extracts boundary splines from splist.
inline std::shared_ptr<PySplineList>
ListExtractBoundaries(const PySplineList& splist,
                      const int nthreads,
                      const bool same_para_dims) {
  const int n_splines = splist.size();
  // to accumulate
  int n_boundaries{};
  // gather offsets for boundary add one last so that we can always find out
  // n_boundary for each spline
  IntVector boundary_offsets{};
  boundary_offsets.reserve(n_splines + 1);

  // in case of same para dim, it'd be equidistance
  int offset{}; // offset counter
  if (same_para_dims) {
    // compute n_boundary
    const int n_boundary = splist[0]->para_dim_ * 2;

    // fill offset

    for (int i{}; i < n_splines; ++i) {
      boundary_offsets.push_back(offset);
      offset += n_boundary;
    }

  } else {
    // for un-equal boundary sizes, we lookup each one of them
    for (int i{}; i < n_splines; ++i) {
      boundary_offsets.push_back(offset);
      offset += splist[i]->para_dim_ * 2;
    }
  }
  // current offset value should equal to total boundary count
  boundary_offsets.push_back(offset);
  n_boundaries = offset;

  // prepare output
  auto out_boundaries = std::make_shared<PySplineList>(n_boundaries);

  // prepare lambda
  auto boundary_extract = [&](int begin, int end) {
    // deref
    auto& ob_deref = *out_boundaries;
    for (int i{begin}; i < end; ++i) {
      // get core
      auto& core = *splist[i]->Core();
      // start of the offset
      auto const& this_offset = boundary_offsets[i];
      // end of the offset
      auto const& next_offset = boundary_offsets[i + 1];
      for (int j{}; j < next_offset - this_offset; ++j) {
        ob_deref[this_offset + j] =
            std::make_shared<PySpline>(core.SplinepyExtractBoundary(j));
      }
    }
  };

  splinepy::utils::NThreadExecution(boundary_extract, n_splines, nthreads);

  return out_boundaries;
}

// computes boundary centers for splines of same para_dim and dim
inline py::array_t<double>
ListBoundaryCenters(const PySplineList& splist,
                    const int nthreads,
                    const bool same_parametric_bounds,
                    const bool check_dims) {
  // dim check first
  if (check_dims) {
    RaiseIfElementsHaveUnequalParaDimOrDim(splist, nthreads);
  }

  // prepare output
  // from here we assume that all the splines have the same para_dim and dim
  const int n_splines = splist.size();
  const int& para_dim = splist[0]->para_dim_;
  const int& dim = splist[0]->dim_;
  const int n_queries = 2 * para_dim;
  const int n_total = n_queries * n_splines;
  py::array_t<double> boundary_centers({n_total, dim});
  double* boundary_centers_ptr =
      static_cast<double*>(boundary_centers.request().ptr);

  // pre-compute boundary centers
  DoubleVector para_bounds;
  double* para_bounds_ptr;
  if (!same_parametric_bounds) {
    para_bounds.resize(n_total * para_dim);
    para_bounds_ptr = para_bounds.data();
    const int stride = n_queries * para_dim;

    auto calc_para_bounds = [&](int begin, int end) {
      for (int i{begin}; i < end; ++i) {
        splinepy::splines::helpers::ScalarTypeBoundaryCenters(
            *splist[i]->Core(),
            &para_bounds_ptr[stride * i]);
      }
    };

    // exe
    splinepy::utils::NThreadExecution(calc_para_bounds, n_splines, nthreads);
  }

  auto calc_boundary_centers_step = [&](int begin, int total_) {
    // each thread needs one query
    DoubleVector queries_vector; /* unused if same_parametric_bounds=true*/
    double* queries;

    // pre compute boundary centers if para bounds are the same
    if (same_parametric_bounds) {
      queries_vector.resize(2 * para_dim * para_dim);
      queries = queries_vector.data();
      splinepy::splines::helpers::ScalarTypeBoundaryCenters(*splist[0]->Core(),
                                                            queries);
    }

    for (int i{begin}; i < total_; i += nthreads) {
      const auto [i_spline, i_query] = std::div(i, n_queries);
      const auto& core = *splist[i_spline]->Core();

      // get ptr start
      if (!same_parametric_bounds) {
        queries = &para_bounds_ptr[i_spline * n_queries * para_dim];
      }

      // eval
      core.SplinepyEvaluate(
          &queries[i_query * para_dim],
          &boundary_centers_ptr[(i_spline * n_queries + i_query) * dim]);
    }
  };

  splinepy::utils::NThreadExecution(calc_boundary_centers_step,
                                    n_total,
                                    nthreads,
                                    splinepy::utils::NThreadQueryType::Step);

  return boundary_centers;
}

inline std::shared_ptr<PySplineList> ListCompose() {}
inline std::shared_ptr<PySplineList> ListCompositionDerivative() {}

/// bind vector of PySpline and add some deprecated cpp functions that maybe
/// nice to have
inline void add_spline_list_pyclass(py::module& m) {

  // use shared_ptr as holder
  auto klasse = py::bind_vector<PySplineList, std::shared_ptr<PySplineList>>(
      m,
      "SplineList");

  klasse
      .def("reserve",
           [](PySplineList& splist, py::ssize_t size) { splist.reserve(size); })
      .def("resize",
           [](PySplineList& splist, py::ssize_t size) { splist.resize(size); })
      .def("shrink_to_fit",
           [](PySplineList& splist) { splist.shrink_to_fit(); })
      .def("swap",
           [](PySplineList& splist_a, PySplineList& splist_b) {
             splist_a.swap(splist_b);
           })
      .def("raise_dim_mismatch",
           &RaiseIfElementsHaveUnequalParaDimOrDim,
           py::arg("nthreads"))
      .def("evaluate",
           &ListEvaluate,
           py::arg("queries"),
           py::arg("nthreads"),
           py::arg("check_dims"))
      .def("sample",
           &ListSample,
           py::arg("resolution"),
           py::arg("nthreads"),
           py::arg("same_parametric_bounds"),
           py::arg("check_dims"))
      .def("extract_boundaries",
           &ListExtractBoundaries,
           py::arg("nthreads"),
           py::arg("same_para_dims"))
      .def("boundary_centers",
           &ListBoundaryCenters,
           py::arg("nthreads"),
           py::arg("same_parametric_bounds"),
           py::arg("check_dims"));
}

} // namespace splinepy::py
