#pragma once

#include <cstdlib>
#include <iostream>
#include <tuple>

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

struct PySplineListNThreadExecutionHelper {
  // thread level info
  int n_common_, n_rem_, spline_start_, spline_end_, query_offset_, n_queries_,
      dim_;
  // spline level info
  int query_start_, query_end_, output_offset_;

  PySplineListNThreadExecutionHelper(const int& begin_id,
                                     const int& end_id,
                                     const int& n_splines,
                                     const int& n_queries,
                                     const int& dim) {
    const int thread_load = end_id - begin_id;

    // number of queries common to all splines
    n_common_ = thread_load / n_splines;

    // number of remaining queries.
    n_rem_ = thread_load % n_splines;

    // beginning spline id of this load
    spline_start_ = begin_id % n_splines;

    // ending spline id of this load
    spline_end_ = (end_id - 1) % n_splines;

    // offset to the very first query for this thread
    query_offset_ = begin_id / n_splines;

    // copy n_queries and dim
    n_queries_ = n_queries;
    dim_ = dim;
  }

  /// computes information required for each spline within the thread.
  void SetSplineId(const int& spline_id) {

    // output ptr offset
    output_offset_ = spline_id * n_queries_ * dim_;

    // init start
    query_start_ = query_offset_;

    // get query start and end
    int n_additional{};
    if (spline_id < spline_start_) {
      n_additional -= 1;
      query_start_ += 1;
    }
    if (n_rem_ != 0 && spline_id <= spline_end_) {
      n_additional += 1;
    }
    query_end_ = query_start_ + n_common_ + n_additional;
  }
};

inline void RaiseIfElementsHaveUnequalParaDimOrDim(const PySplineList& splist,
                                                   const int nthreads) {
  // use first spline as guide line
  const int& para_dim = splist[0]->para_dim_;
  const int& dim = splist[0]->dim_;

  auto check_dims = [&](int begin, int end) {
    for (int i{begin}; i < end; ++i) {
      CheckParaDimAndDim(*splist[i], para_dim, dim, true);
    }
  };

  splinepy::utils::NThreadExecution(check_dims,
                                    static_cast<int>(splist.size()),
                                    nthreads);
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
  auto evaluate = [&](int begin, int end) {
    auto thread_helper = PySplineListNThreadExecutionHelper(begin,
                                                            end,
                                                            n_splines,
                                                            n_queries,
                                                            dim);

    // loop splines
    for (int i{}; i < n_splines; ++i) {
      const auto& spl = *splist[i];
      const auto& core = *spl.Core();

      // compute start and end of query, and output ptr offset.
      thread_helper.SetSplineId(i);

      // queries for splines
      for (int j{thread_helper.query_start_}; j < thread_helper.query_end_;
           ++j) {
        core.SplinepyEvaluate(
            &queries_ptr[j * para_dim],
            &evaluated_ptr[thread_helper.output_offset_ + (j * dim)]);
      }
    }
  };

  // exe
  splinepy::utils::NThreadExecution(evaluate, n_total, nthreads);

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

  // queries - will only be filled if same_parametric_bounds=true
  DoubleVector queries_vector;
  double* queries;

  // create variable for lambda - needs different ones based on
  // same_parametric_bounds
  std::function<void(int, int)> thread_func;

  // GridPoints for same_parametric_bounds=true, else, won't be used
  splinepy::utils::DefaultInitializationVector<
      splinepy::utils::CStyleArrayPointerGridPoints>
      grid_points;

  // if you know all the queries have same parametric bounds
  // you don't need to re-compute queries
  if (same_parametric_bounds) {
    // get para bounds
    DoubleVector para_bounds_vector(2 * para_dim);
    double* para_bounds = para_bounds_vector.data();
    first_spline.Core()->SplinepyParametricBounds(para_bounds);

    splinepy::utils::CStyleArrayPointerGridPoints gp_generator(para_dim,
                                                               para_bounds,
                                                               resolutions);
    // assign queries
    queries_vector.resize(n_queries * para_dim);
    queries = queries_vector.data();

    gp_generator.Fill(queries);

    // create lambda for nthread exe
    thread_func = [&](int begin, int end) {
      auto thread_helper = PySplineListNThreadExecutionHelper(begin,
                                                              end,
                                                              n_splines,
                                                              n_queries,
                                                              dim);

      // loop splines
      for (int i{}; i < n_splines; ++i) {
        const auto& core = *splist[i]->Core();

        // compute start and end of query, and output ptr offset.
        thread_helper.SetSplineId(i);

        // queries for splines
        for (int j{thread_helper.query_start_}; j < thread_helper.query_end_;
             ++j) {
          core.SplinepyEvaluate(
              &queries[j * para_dim],
              &sampled_ptr[thread_helper.output_offset_ + (j * dim)]);
        }
      }
    };

  } else {
    // we assume each spline has different parametric bounds

    // resize grid points to have same size as input spline list
    grid_points.resize(n_splines);

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

    // pre compute entries
    splinepy::utils::NThreadExecution(create_grid_points, n_splines, nthreads);

    // similar to the one with same_parametric_bounds, except it computes query
    // on the fly
    thread_func = [&](int begin, int end) {
      auto thread_helper = PySplineListNThreadExecutionHelper(begin,
                                                              end,
                                                              n_splines,
                                                              n_queries,
                                                              dim);
      // each thread needs just one query array
      DoubleVector thread_query_vector(para_dim);
      double* thread_query = thread_query_vector.data();
      // loop splines
      for (int i{}; i < n_splines; ++i) {
        // get spline core and grid point helper
        const auto& core = *splist[i]->Core();
        const auto& gp_helper = grid_points[i];

        // compute start and end of query, and output ptr offset.
        thread_helper.SetSplineId(i);

        // queries for splines
        for (int j{thread_helper.query_start_}; j < thread_helper.query_end_;
             ++j) {
          gp_helper.IdToGridPoint(j, thread_query);

          core.SplinepyEvaluate(
              thread_query,
              &sampled_ptr[thread_helper.output_offset_ + (j * dim)]);
        }
      }
    };
  }

  // nthread exe
  splinepy::utils::NThreadExecution(thread_func, n_total, nthreads);

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
  const int n_splines = splist.size();
  const int& para_dim = splist[0]->para_dim_;
  const int& dim = splist[0]->dim_;
  const int n_queries = 2 * para_dim;
  const int n_total = n_queries * n_splines;
  py::array_t<double> boundary_centers({n_total, dim});
  double* boundary_centers_ptr =
      static_cast<double*>(boundary_centers.request().ptr);

  // we assume from here that all the splines have the same para_dim and dim
  auto calc_boundary_centers = [&](int begin, int end) {
    // each thread needs one query
    DoubleVector queries_vector(2 * para_dim * para_dim);
    double* queries = queries_vector.data();

    // pre compute boundary centers if para bounds are the same
    if (same_parametric_bounds) {
      splinepy::splines::helpers::ScalarTypeBoundaryCenters(*splist[0]->Core(),
                                                            queries);
    }

    // prepare thread helper
    auto thread_helper = PySplineListNThreadExecutionHelper(begin,
                                                            end,
                                                            n_splines,
                                                            n_queries,
                                                            dim);
    // loop splines
    for (int i{}; i < n_splines; ++i) {
      // get core spline (SplinepyBase)
      const auto& core = *splist[i]->Core();

      // compute start and end of query, and output ptr offset.
      thread_helper.SetSplineId(i);

      // in case para_bounds are assumed to be different, compute query
      if (!same_parametric_bounds) {
        splinepy::splines::helpers::ScalarTypeBoundaryCenters(core, queries);
      }

      // queries for splines
      for (int j{thread_helper.query_start_}; j < thread_helper.query_end_;
           ++j) {
        core.SplinepyEvaluate(
            &queries[j * para_dim],
            &boundary_centers_ptr[thread_helper.output_offset_ + (j * dim)]);
      }
    }
  };

  splinepy::utils::NThreadExecution(calc_boundary_centers, n_total, nthreads);

  return boundary_centers;
}

/// bind vector of PySpline and add some deprecated cpp functions that maybe
/// nice to have
inline void add_spline_list_pyclass(py::module& m) {

  // use shared_ptr as holder
  py::bind_vector<PySplineList, std::shared_ptr<PySplineList>>(m, "SplineList");

  m.def("list_reserve",
        [](PySplineList& splist, py::ssize_t size) { splist.reserve(size); });
  m.def("list_resize",
        [](PySplineList& splist, py::ssize_t size) { splist.resize(size); });
  m.def("list_shrink_to_fit",
        [](PySplineList& splist) { splist.shrink_to_fit(); });
  m.def("list_swap", [](PySplineList& splist_a, PySplineList& splist_b) {
    splist_a.swap(splist_b);
  });
  m.def("list_evaluate",
        &ListEvaluate,
        py::arg("spline_list"),
        py::arg("queries"),
        py::arg("nthreads"),
        py::arg("check_dims"));
  m.def("list_sample",
        &ListSample,
        py::arg("spline_list"),
        py::arg("resolution"),
        py::arg("nthreads"),
        py::arg("same_parametric_bounds"),
        py::arg("check_dims"));
  m.def("list_extract_boundaries",
        &ListExtractBoundaries,
        py::arg("spline_list"),
        py::arg("n_threads"),
        py::arg("same_para_dims"));
  m.def("list_boundary_centers",
        &ListBoundaryCenters,
        py::arg("spline_list"),
        py::arg("nthreads"),
        py::arg("same_parametric_bounds"),
        py::arg("check_dims"));
}

} // namespace splinepy::py
