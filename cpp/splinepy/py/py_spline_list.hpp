#pragma once

#include <cstdlib>
#include <iostream>
#include <tuple>

// pybind
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

#include <splinepy/py/py_spline.hpp>
#include <splinepy/py/py_spline_extensions.hpp>
#include <splinepy/utils/nthreads.hpp>

// PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<splinepy::py::PySpline>>);

namespace splinepy::py {

namespace py = pybind11;

using PySplineList = std::vector<std::shared_ptr<PySpline>>;

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
inline py::array_t<double> EvaluateList(const PySplineList& splist,
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

inline py::array_t<double> DerivativeList(const PySplineList& splist,
                                          const py::array_t<double>& queries,
                                          const py::array_t<int>& orders,
                                          const int nthreads) {}

/// Samples equal resoltions for each para dim
inline py::array_t<double> SampleList(const PySplineList& splist,
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
  int* resolutions = new int[para_dim];
  for (int i{}; i < para_dim; ++i) {
    resolutions[i] = resolution;
    n_queries *= resolution;
  }

  // prepare input /  output
  const int n_total = n_splines * n_queries;
  py::array_t<double> sampled({n_total, dim});
  double* sampled_ptr = static_cast<double*>(sampled.request().ptr);

  // queries - will only be filled if same_parametric_bounds=true
  double* queries;

  // create variable for lambda - needs different ones based on
  // same_parametric_bounds
  std::function<void(int, int)> thread_func;

  // GridPoints for same_parametric_bounds=true, this will have size=1.
  std::vector<splinepy::utils::CStyleArrayPointerGridPoints> grid_points;

  // if you know all the queries have same parametric bounds
  // you don't need to re-compute queries
  if (same_parametric_bounds) {
    // get para bounds
    double* para_bounds = new double[2 * para_dim];
    first_spline.Core()->SplinepyParametricBounds(para_bounds);

    // create grid points helper
    grid_points.emplace_back(
        splinepy::utils::CStyleArrayPointerGridPoints(para_dim,
                                                      para_bounds,
                                                      resolutions));
    const auto& gp_generator = grid_points[0];

    // assign queries
    queries = new double[n_queries * para_dim];

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

    // release
    delete[] para_bounds;

  } else {
    // we assume each spline has different parametric bounds

    // resize grid points to have same size as input spline list
    grid_points.resize(n_splines);

    // create grid_points
    auto create_grid_points = [&](int begin, int end) {
      double* para_bounds = new double[2 * para_dim];

      for (int i{begin}; i < end; ++i) {
        // get para_bounds
        splist[i]->Core()->SplinepyParametricBounds(para_bounds);
        // setup grid points helper
        grid_points[i].SetUp(para_dim, para_bounds, resolutions);
      }
      delete[] para_bounds;
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
      double* q_ptr = new double[para_dim];
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
          gp_helper.IdToGridPoint(j, q_ptr);

          core.SplinepyEvaluate(
              q_ptr,
              &sampled_ptr[thread_helper.output_offset_ + (j * dim)]);
        }
      }
      delete[] q_ptr;
    };
  }

  // nthread exe
  splinepy::utils::NThreadExecution(thread_func, n_total, nthreads);

  // release
  delete[] resolutions;
  if (same_parametric_bounds) {
    delete[] queries;
  }

  return sampled;
}

inline std::shared_ptr<PySplineList>
ExtractBoundarySplines(const PySplineList& splist, const int nthreads) {}

inline py::array_t<double> BoundaryCenters(const PySplineList& splist,
                                           const int nthreads) {}

/// bind vector of PySpline and add some deprecated cpp functions that maybe
/// nice to have
inline void add_spline_list_pyclass(py::module& m) {

  // use shared_ptr as holder
  py::bind_vector<PySplineList, std::shared_ptr<PySplineList>>(m, "SplineList");

  m.def("reserve_list",
        [](PySplineList& splist, py::ssize_t size) { splist.reserve(size); });
  m.def("resize_list",
        [](PySplineList& splist, py::ssize_t size) { splist.resize(size); });
  m.def("shrink_to_fit_list",
        [](PySplineList& splist) { splist.shrink_to_fit(); });
  m.def("swap_list", [](PySplineList& splist_a, PySplineList& splist_b) {
    splist_a.swap(splist_b);
  });
  m.def("evaluate_list",
        &EvaluateList,
        py::arg("spline_list"),
        py::arg("queries"),
        py::arg("nthreads"),
        py::arg("check_dims"));
  m.def("sample_list",
        &SampleList,
        py::arg("spline_list"),
        py::arg("resolution"),
        py::arg("nthreads"),
        py::arg("same_parametric_bounds"),
        py::arg("check_dims"));
}

} // namespace splinepy::py
