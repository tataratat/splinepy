#pragma once

#include <cstdlib>
#include <string>

// pybind
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

#include <splinepy/py/py_spline.hpp>
#include <splinepy/py/py_spline_extensions.hpp>
#include <splinepy/splines/helpers/scalar_type_wrapper.hpp>
#include <splinepy/splines/splinepy_base.hpp>
#include <splinepy/utils/default_initialization_allocator.hpp>
#include <splinepy/utils/nthreads.hpp>

namespace splinepy::py {

namespace py = pybind11;

// alias
using CoreSplineVector =
    std::vector<std::shared_ptr<splinepy::splines::SplinepyBase>>;
using IntVector = splinepy::utils::DefaultInitializationVector<int>;
using IntVectorVector = splinepy::utils::DefaultInitializationVector<IntVector>;
using DoubleVector = splinepy::utils::DefaultInitializationVector<double>;

/// Internal use only
/// Extract CoreSpline_s from list of PySplines
inline std::vector<PySpline::CoreSpline_>
ListOfPySplinesToVectorOfCoreSplines(py::list pysplines,
                                     const int nthreads = 1) {
  // prepare return obj
  const int n_splines = static_cast<int>(pysplines.size());
  std::vector<PySpline::CoreSpline_> core_splines(n_splines);

  auto to_core = [&](int begin, int end) {
    for (int i{begin}; i < end; ++i) {
      core_splines[i] =
          pysplines[i].template cast<std::shared_ptr<PySpline>>()->Core();
    }
  };
  splinepy::utils::NThreadExecution(to_core, n_splines, nthreads);

  return core_splines;
}

/// TODO: move ListOfPySplinesToVectorOfCoreSplines here and rename
inline py::list CoreSplineVectorToPySplineList(CoreSplineVector& splist,
                                               int nthreads) {
  // prepare return obj
  const int n_splines = static_cast<int>(splist.size());
  py::list pyspline_list(n_splines);

  auto to_pyspline = [&](int begin, int end) {
    for (int i{begin}; i < end; ++i) {
      pyspline_list[i] = std::make_shared<PySpline>(splist[i]);
    }
  };

  // multi thread execution causes segfault.
  // until we find a better solution, do single exe
  nthreads = 1;

  splinepy::utils::NThreadExecution(to_pyspline, n_splines, nthreads);

  return pyspline_list;
}

/// @brief raises if elements' para_dim and dims aren't equal
/// @param splist
/// @param nthreads
inline void
RaiseIfElementsHaveUnequalParaDimOrDim(const CoreSplineVector& splist,
                                       const int nthreads) {
  // use first spline as guide line
  const int para_dim = splist[0]->SplinepyParaDim();
  const int dim = splist[0]->SplinepyDim();

  std::vector<IntVector> mismatches(nthreads);

  auto check_dims_step = [&](int begin, int total_) {
    for (int i{begin}; i < total_; i += nthreads) {
      if (!CheckCoreParaDimAndDim(splist[i], para_dim, dim, false)) {
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

inline void ListRaiseIfElementsHaveUnequalParaDimOrDim(py::list& splist,
                                                       const int nthreads) {
  const auto core_vector =
      ListOfPySplinesToVectorOfCoreSplines(splist, nthreads);

  RaiseIfElementsHaveUnequalParaDimOrDim(core_vector, nthreads);
}

/// evaluates splines of same para_dim and dim.
inline py::array_t<double>
CoreSplineVectorEvaluate(const CoreSplineVector& splist,
                         const py::array_t<double>& queries,
                         const int nthreads,
                         const bool check_dims) {
  // use first spline as dimension guide line
  const int para_dim = splist[0]->SplinepyParaDim();
  const int dim = splist[0]->SplinepyDim();

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
      splist[i_spline]->SplinepyEvaluate(
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

inline py::array_t<double> ListEvaluate(py::list& splist,
                                        const py::array_t<double>& queries,
                                        const int nthreads,
                                        const bool check_dims) {
  const auto core_vector =
      ListOfPySplinesToVectorOfCoreSplines(splist, nthreads);

  return CoreSplineVectorEvaluate(core_vector, queries, nthreads, check_dims);
}

inline py::array_t<double>
CoreSplineVectorDerivative(const CoreSplineVector& splist,
                           const py::array_t<double>& queries,
                           const py::array_t<int>& orders,
                           const int nthreads) {
  return py::array_t<double>();
}

inline py::array_t<double> ListDerivative(py::list& splist,
                                          const py::array_t<double>& queries,
                                          const py::array_t<int>& orders,
                                          const int nthreads) {
  const auto core_vector =
      ListOfPySplinesToVectorOfCoreSplines(splist, nthreads);

  return CoreSplineVectorDerivative(core_vector, queries, orders, nthreads);
}

/// Samples equal resoltions for each para dim
inline py::array_t<double>
CoreSplineVectorSample(const CoreSplineVector& splist,
                       const int resolution,
                       const int nthreads,
                       const bool same_parametric_bounds,
                       const bool check_dims) {
  // use first spline as dimension guide line
  const auto& first_spline = *splist[0];
  const int para_dim = first_spline.SplinepyParaDim();
  const int dim = first_spline.SplinepyDim();

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
    first_spline.SplinepyParametricBounds(para_bounds);

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
        splist[i_spline]->SplinepyEvaluate(
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
        splist[i]->SplinepyParametricBounds(para_bounds);
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
        splist[i_spline]->SplinepyEvaluate(
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

/// Samples equal resoltions for each para dim
inline py::array_t<double> ListSample(py::list& splist,
                                      const int resolution,
                                      const int nthreads,
                                      const bool same_parametric_bounds,
                                      const bool check_dims) {
  const auto core_vector =
      ListOfPySplinesToVectorOfCoreSplines(splist, nthreads);

  return CoreSplineVectorSample(core_vector,
                                resolution,
                                nthreads,
                                same_parametric_bounds,
                                check_dims);
}

// extracts boundary splines from splist.
inline CoreSplineVector
CoreSplineVectorExtractBoundaries(const CoreSplineVector& splist,
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
    const int n_boundary = splist[0]->SplinepyParaDim() * 2;

    // fill offset

    for (int i{}; i < n_splines; ++i) {
      boundary_offsets.push_back(offset);
      offset += n_boundary;
    }

  } else {
    // for un-equal boundary sizes, we lookup each one of them
    for (int i{}; i < n_splines; ++i) {
      boundary_offsets.push_back(offset);
      offset += splist[i]->SplinepyParaDim() * 2;
    }
  }
  // current offset value should equal to total boundary count
  boundary_offsets.push_back(offset);
  n_boundaries = offset;

  // prepare output
  CoreSplineVector out_boundaries(n_boundaries);

  // prepare lambda
  auto boundary_extract = [&](int begin, int end) {
    for (int i{begin}; i < end; ++i) {
      // start of the offset
      const auto& this_offset = boundary_offsets[i];
      // end of the offset
      const auto& next_offset = boundary_offsets[i + 1];
      // get spline
      auto& spline_i = *splist[i];
      for (int j{}; j < next_offset - this_offset; ++j) {
        out_boundaries[this_offset + j] = spline_i.SplinepyExtractBoundary(j);
      }
    }
  };

  splinepy::utils::NThreadExecution(boundary_extract, n_splines, nthreads);

  return out_boundaries;
}

inline py::list ListExtractBoundaries(py::list& splist,
                                      const int nthreads,
                                      const bool same_para_dims) {
  const auto core_vector =
      ListOfPySplinesToVectorOfCoreSplines(splist, nthreads);
  auto boundaries =
      CoreSplineVectorExtractBoundaries(core_vector, nthreads, same_para_dims);

  return CoreSplineVectorToPySplineList(boundaries, nthreads);
}

// computes boundary centers for splines of same para_dim and dim
inline py::array_t<double>
CoreSplineVectorBoundaryCenters(const CoreSplineVector& splist,
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
  const int para_dim = splist[0]->SplinepyParaDim();
  const int dim = splist[0]->SplinepyDim();
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
            *splist[i],
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
      splinepy::splines::helpers::ScalarTypeBoundaryCenters(*splist[0],
                                                            queries);
    }

    for (int i{begin}; i < total_; i += nthreads) {
      const auto [i_spline, i_query] = std::div(i, n_queries);

      // get ptr start
      if (!same_parametric_bounds) {
        queries = &para_bounds_ptr[i_spline * n_queries * para_dim];
      }

      // eval
      splist[i_spline]->SplinepyEvaluate(
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

inline py::array_t<double>
ListBoundaryCenters(py::list& splist,
                    const int nthreads,
                    const bool same_parametric_bounds,
                    const bool check_dims) {
  const auto core_vector =
      ListOfPySplinesToVectorOfCoreSplines(splist, nthreads);
  return CoreSplineVectorBoundaryCenters(core_vector,
                                         nthreads,
                                         same_parametric_bounds,
                                         check_dims);
}

/// @brief NThread Compose. Invalid inputs will raise runtime_error on the fly
/// @param outer_splines
/// @param inner_splines
/// @param cartesian_product If true, composes each inner splines to all outer
/// splines. Else, outer and inner splines must have same len.
/// @param nthreads
/// @return
inline CoreSplineVector
CoreSplineVectorCompose(const CoreSplineVector& outer_splines,
                        const CoreSplineVector& inner_splines,
                        const bool cartesian_product,
                        const int nthreads) {

  // get size
  const int n_outer = outer_splines.size();
  const int n_inner = inner_splines.size();

  // create output
  CoreSplineVector composed_splines;

  // check if size matchs
  if (!cartesian_product) {
    if (n_outer != n_inner) {
      splinepy::utils::PrintAndThrowError(
          "Length mismatch of outer_splines (",
          n_outer,
          ") and inner_splines (",
          n_inner,
          "). To compose each inner splines to all outer splines, please set "
          "cartesian_product=True.");
    }

    // acquisition of output space
    composed_splines.resize(n_outer);

    auto compose = [&](int begin, int end) {
      for (int i{begin}; i < end; ++i) {
        composed_splines[i] =
            outer_splines[i]->SplinepyCompose(inner_splines[i]);
      }
    };

    splinepy::utils::NThreadExecution(compose, n_outer, nthreads);

  } else {
    const int n_total = n_outer * n_inner;
    composed_splines.resize(n_total);

    auto compose_step = [&](int begin, int total_) {
      for (int i{begin}; i < total_; ++i) {
        const auto [i_outer, i_inner] = std::div(i, n_inner);
        composed_splines[i] =
            outer_splines[i_outer]->SplinepyCompose(inner_splines[i_inner]);
      }
    };
    splinepy::utils::NThreadExecution(compose_step,
                                      n_total,
                                      nthreads,
                                      splinepy::utils::NThreadQueryType::Step);
  }

  return composed_splines;
}

inline py::list ListCompose(py::list& outer_splines,
                            py::list& inner_splines,
                            const bool cartesian_product,
                            const int nthreads) {
  const auto core_outer =
      ListOfPySplinesToVectorOfCoreSplines(outer_splines, nthreads);
  const auto core_inner =
      ListOfPySplinesToVectorOfCoreSplines(inner_splines, nthreads);

  auto composed_cores = CoreSplineVectorCompose(core_outer,
                                                core_inner,
                                                cartesian_product,
                                                nthreads);

  return CoreSplineVectorToPySplineList(composed_cores, nthreads);
}

/// @brief NThread composition derivative. same query options as ListCompose()
/// @param outer_splines
/// @param inner_splines
/// @param inner_derivative
/// @param cartesian_product
/// @param nthreads
/// @return
inline CoreSplineVector
CoreSplineVectorCompositionDerivative(const CoreSplineVector& outer_splines,
                                      const CoreSplineVector& inner_splines,
                                      const CoreSplineVector& inner_derivatives,
                                      const bool cartesian_product,
                                      const int nthreads) {
  // get size
  const int n_outer = outer_splines.size();
  const int n_inner = inner_splines.size();
  const int n_inner_der = inner_derivatives.size();

  // number of queries equals to para_dim
  const int& n_queries = inner_splines[0]->SplinepyParaDim();

  // create output
  // create derivative spline container in case of cartesian product
  CoreSplineVector composition_derivatives, outer_spline_derivatives,
      inner_derivatives_single_dims;

  // check if size matchs
  // 1. inner and inner_der
  if (n_inner != n_inner_der) {
    splinepy::utils::PrintAndThrowError(
        "Length mismatch - inner_splines (",
        n_inner,
        ") and inner_derivatives (",
        n_inner_der,
        ") should have same length. To compose each inner splines to all outer "
        "splines, please set cartesian_product=True.");
  }
  // 2. outer and inner incase of not - cartesian product.
  if (!cartesian_product) {
    if (n_outer != n_inner) {
      splinepy::utils::PrintAndThrowError(
          "Length mismatch of outer_splines (",
          n_outer,
          ") and inner_splines (",
          n_inner,
          "). To compose each inner splines to all outer splines, please set "
          "cartesian_product=True.");
    }

    // acquisition of output space
    composition_derivatives.resize(n_outer);

    // exe
    auto calc_composition_derivatives = [&](int begin, int end) {
      for (int i{begin}; i < end; ++i) {
        composition_derivatives[i] =
            outer_splines[i]->SplinepyCompositionDerivative(
                inner_splines[i],
                inner_derivatives[i]);
      }
    };

    splinepy::utils::NThreadExecution(calc_composition_derivatives,
                                      n_outer,
                                      nthreads);

  } else {
    const int n_total = n_outer * n_inner;
    composition_derivatives.resize(n_total);

    const int n_outer_precompute = n_outer * n_queries;
    const int n_inner_precompute = n_inner * n_queries;

    // resize vectors to hold der splines of outer splines and der splines
    // their para_dim and dim should match. Otherwise, it will raise.
    outer_spline_derivatives.resize(n_outer_precompute);
    inner_derivatives_single_dims.resize(n_inner_precompute);

    // make query arrays
    IntVectorVector order_queries(n_queries);
    for (int i{}; i < n_queries; ++i) {
      auto& order_query = order_queries[i];
      order_query.resize(n_queries, 0);
      order_query[i] = 1;
    }

    // let's fill der splines and extracted splines
    auto calc_der_splines_and_extract_dim = [&](int begin, int total_) {
      for (int i{begin}; i < total_; i += nthreads) {
        const auto [i_spline, i_query] = std::div(i, n_queries);

        const int offset = i_spline * n_queries;
        // fill if outer is still in range
        if (i_spline < n_outer) {
          outer_spline_derivatives[offset + i_query] =
              outer_splines[i_spline]->SplinepyDerivativeSpline(
                  order_queries[i_query].data());
        }

        if (i_spline < n_inner) {
          inner_derivatives_single_dims[offset + i_query] =
              inner_derivatives[i_spline]->SplinepyExtractDim(i_query);
        }
      }
    };

    // precompute
    const int precompute_total =
        std::max(n_outer_precompute, n_inner_precompute);
    splinepy::utils::NThreadExecution(calc_der_splines_and_extract_dim,
                                      precompute_total,
                                      nthreads,
                                      splinepy::utils::NThreadQueryType::Step);
    // now, composition der
    auto calc_composition_derivatives_step = [&](int begin, int total_) {
      for (int i{begin}; i < total_; i += nthreads) {
        const auto [i_outer, i_inner] = std::div(i, n_inner);
        // frequently used core
        const auto& inner_core = inner_splines[i_inner];

        // this one needs a loop
        // create
        auto this_comp_der =
            outer_spline_derivatives[i_outer * n_queries]
                ->SplinepyCompose(inner_core)
                ->SplinepyMultiply(
                    inner_derivatives_single_dims[i_inner * n_queries]);

        // add
        for (int j{1}; j < n_queries; ++j) {
          this_comp_der = this_comp_der->SplinepyAdd(
              outer_spline_derivatives[i_outer * n_queries + j]
                  ->SplinepyCompose(inner_core)
                  ->SplinepyMultiply(
                      inner_derivatives_single_dims[i_inner * n_queries + j]));
        }
        // now fill
        composition_derivatives[i] = this_comp_der;
      }
    };

    splinepy::utils::NThreadExecution(calc_composition_derivatives_step,
                                      n_total,
                                      nthreads,
                                      splinepy::utils::NThreadQueryType::Step);
  }

  return composition_derivatives;
}

inline py::list ListCompositionDerivative(py::list& outer_splines,
                                          py::list& inner_splines,
                                          py::list& inner_derivatives,
                                          const bool cartesian_product,
                                          const int nthreads) {
  const auto core_outer =
      ListOfPySplinesToVectorOfCoreSplines(outer_splines, nthreads);
  const auto core_inner =
      ListOfPySplinesToVectorOfCoreSplines(inner_splines, nthreads);
  const auto core_inner_derivatives =
      ListOfPySplinesToVectorOfCoreSplines(inner_derivatives, nthreads);

  auto core_com_der =
      CoreSplineVectorCompositionDerivative(core_outer,
                                            core_inner,
                                            core_inner_derivatives,
                                            cartesian_product,
                                            nthreads);

  return CoreSplineVectorToPySplineList(core_com_der, nthreads);
}

/// bind vector of PySpline and add some deprecated cpp functions that maybe
/// nice to have
inline void add_spline_list_pyclass(py::module& m) {

  auto list_module = m.def_submodule(
      "lists",
      "Module for list based queries with multithreading capabilities. Please "
      "make sure splines don't have any inplace modifications.");

  list_module
      .def("raise_dim_mismatch",
           &ListRaiseIfElementsHaveUnequalParaDimOrDim,
           py::arg("spline_list"),
           py::arg("nthreads"))
      .def("evaluate",
           &ListEvaluate,
           py::arg("spline_list"),
           py::arg("queries"),
           py::arg("nthreads"),
           py::arg("check_dims"))
      .def("sample",
           &ListSample,
           py::arg("spline_list"),
           py::arg("resolution"),
           py::arg("nthreads"),
           py::arg("same_parametric_bounds"),
           py::arg("check_dims"))
      .def("extract_boundaries",
           &ListExtractBoundaries,
           py::arg("spline_list"),
           py::arg("nthreads"),
           py::arg("same_para_dims"))
      .def("boundary_centers",
           &ListBoundaryCenters,
           py::arg("spline_list"),
           py::arg("nthreads"),
           py::arg("same_parametric_bounds"),
           py::arg("check_dims"))
      .def("compose",
           &ListCompose,
           py::arg("outer_splines"),
           py::arg("inner_splines"),
           py::arg("cartesian_product"),
           py::arg("nthreads"))
      .def("composition_derivative",
           &ListCompositionDerivative,
           py::arg("outer_splines"),
           py::arg("inner_splines"),
           py::arg("inner_derivatives"),
           py::arg("cartesian_product"),
           py::arg("nthreads"));
}

} // namespace splinepy::py
