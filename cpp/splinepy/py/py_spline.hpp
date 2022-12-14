#pragma once

#include <algorithm>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

// pybind11
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

// first four are required for Create* implmentations
#include <splinepy/splines/bezier.hpp>
#include <splinepy/splines/bspline.hpp>
#include <splinepy/splines/nurbs.hpp>
#include <splinepy/splines/rational_bezier.hpp>
#include <splinepy/splines/splinepy_base.hpp>
#include <splinepy/utils/grid_points.hpp>
#include <splinepy/utils/print.hpp>

namespace splinepy::py {

namespace py = pybind11;

using namespace splinelib::sources;

template<typename ValueType>
static bool CheckPyArrayShape(const py::array_t<ValueType> arr,
                              const std::vector<int>& shape,
                              const bool throw_ = true) {
  const std::size_t expected_dim = shape.size();
  if (expected_dim != static_cast<std::size_t>(arr.ndim())) {
    if (!throw_)
      return false;
    splinepy::utils::PrintAndThrowError("Array dim mismatch.",
                                        "Expected -",
                                        expected_dim,
                                        "Given -",
                                        arr.ndim());
  }
  const py::ssize_t* arrshape = arr.shape();
  for (std::size_t i{}; i < expected_dim; ++i) {
    const int& shape_i = shape[i];
    if (shape_i < 0) {
      continue;
    } else {
      if (shape_i != arrshape[i]) {
        if (!throw_)
          return false;
        splinepy::utils::PrintAndThrowError("Array shape mismatch",
                                            "in dimension [",
                                            i,
                                            "].",
                                            "Expected -",
                                            shape_i,
                                            "Given -",
                                            arrshape[i]);
      }
    }
  }
  return true;
}

template<typename ValueType>
static bool CheckPyArraySize(const py::array_t<ValueType> arr,
                             const int size,
                             const bool throw_ = true) {
  if (size != arr.size()) {
    if (!throw_)
      return false;
    splinepy::utils::PrintAndThrowError("array size mismatch",
                                        "expected -",
                                        size,
                                        "given -",
                                        arr.size());
  }
  return true;
}

/// True interface to python
class PySpline {
public:
  using CoreSpline_ = typename std::shared_ptr<splinepy::splines::SplinepyBase>;

  CoreSpline_ c_spline_ = nullptr;

  // store very frequently used values - from python this will be readonly
  int para_dim_ = -1;
  int dim_ = -1;

  // place to store and access from both python and cpp side.
  py::dict data_;

  // ctor
  PySpline() = default;
  PySpline(PySpline&&) = default;
  PySpline(const py::kwargs& kwargs) { NewCore(kwargs); }
  PySpline(const CoreSpline_& another_core) : c_spline_(another_core) {
    para_dim_ = c_spline_->SplinepyParaDim();
    dim_ = c_spline_->SplinepyDim();
  }
  PySpline(PySpline& another_py_spline)
      : c_spline_(another_py_spline.Core()),
        para_dim_(another_py_spline.para_dim_),
        dim_(another_py_spline.dim_) {
    // nichts
  }

  /// Creates a corresponding spline based on kwargs
  /// similar to previous update_c()
  /// Runs sanity checks on inputs
  void NewCore(const py::kwargs& kwargs) {
    // parse kwargs
    int* degrees_ptr = nullptr;
    std::vector<std::vector<double>> knot_vectors;
    std::vector<std::vector<double>>* knot_vectors_ptr = nullptr;
    double* control_points_ptr = nullptr;
    double* weights_ptr = nullptr;
    int para_dim, dim, /* for checks -> */ ncps, required_ncps{1};

    // get degrees and set para_dim
    auto d_array = py::cast<py::array_t<int>>(kwargs["degrees"]);
    degrees_ptr = static_cast<int*>(d_array.request().ptr);
    para_dim = d_array.size();
    knot_vectors.reserve(para_dim);

    // get cps and set dim. set ncps for checks later
    auto cp_array = py::cast<py::array_t<double>>(kwargs["control_points"]);
    control_points_ptr = static_cast<double*>(cp_array.request().ptr);
    dim = cp_array.shape(1);
    ncps = cp_array.shape(0);

    // maybe, get knot_vectors
    if (kwargs.contains("knot_vectors")) {
      // check list
      int kv_dim{0};
      double prev_knot, this_knot;

      // loop over list of list/arrays
      for (py::handle kv : kwargs["knot_vectors"]) {
        std::vector<double> knot_vector;
        // cast to array_t, as python side will try to use tracked array anyways
        auto kv_array = py::cast<py::array_t<double>>(kv);
        knot_vector.reserve(kv_array.size());

        int nknots{0};
        prev_knot = -1.;
        for (py::handle k : kv) {
          this_knot = k.cast<double>();

          // can't be negative
          if (this_knot < 0) {
            splinepy::utils::PrintAndThrowError("Parametric dimension (",
                                                kv_dim,
                                                ")",
                                                "includes negative knot.");
          }
          // must be increasing
          if (prev_knot - this_knot > 0) {
            splinepy::utils::PrintAndThrowError(
                prev_knot,
                this_knot,
                "Knots of parametric dimension (",
                kv_dim,
                ")",
                "are not in increasing order.");
          }

          // lgtm, add!
          knot_vector.push_back(this_knot);

          // prepare next
          prev_knot = this_knot;
          ++nknots;
        }

        // minimal number of knots check
        int n_minimal_knots = 2 * (degrees_ptr[kv_dim] + 1);
        if (nknots < n_minimal_knots) {
          splinepy::utils::PrintAndThrowError(
              "Not enough knots in parametric dimension (",
              kv_dim,
              ").",
              "At least",
              n_minimal_knots,
              "knots are expected, but",
              nknots,
              "were given.");
        }
        // multiply expected control mesh resolution
        required_ncps *= nknots - degrees_ptr[kv_dim] - 1;

        // lgtm, add!
        knot_vectors.push_back(std::move(knot_vector));

        // prepare next
        ++kv_dim;
      }

      // check number of control points
      if (ncps != required_ncps) {
        splinepy::utils::PrintAndThrowError("Invalid number of control points.",
                                            required_ncps,
                                            "exepcted, but",
                                            ncps,
                                            "were given.");
      }
      // last, dim check
      if (para_dim != kv_dim) {
        splinepy::utils::PrintAndThrowError(
            "Dimension mis-match between `degrees` (",
            para_dim,
            ") and `knot_vectors` (",
            kv_dim,
            ").");
      }

      // lgtm, set!
      knot_vectors_ptr = &knot_vectors;
    } else {
      // Bezier here. get required ncps
      required_ncps = std::accumulate(degrees_ptr,
                                      degrees_ptr + para_dim,
                                      1,
                                      [](const int& cummulated, const int& d) {
                                        return cummulated * (d + 1);
                                      });
    }

    // check number of control points
    if (ncps != required_ncps) {
      splinepy::utils::PrintAndThrowError("Invalid number of control points.",
                                          required_ncps,
                                          "exepcted, but",
                                          ncps,
                                          "were given.");
    }

    // maybe, get weights
    if (kwargs.contains("weights")) {
      auto w_array = py::cast<py::array_t<double>>(kwargs["weights"]);
      weights_ptr = static_cast<double*>(w_array.request().ptr);

      // check if size matches with ncps
      if (w_array.size() != ncps) {
        splinepy::utils::PrintAndThrowError(
            "Number of weights (",
            w_array.size(),
            ") does not match number of control points (",
            ncps,
            ").");
      }
    }

    // new assign
    c_spline_ =
        splinepy::splines::SplinepyBase::SplinepyCreate(para_dim,
                                                        dim,
                                                        degrees_ptr,
                                                        knot_vectors_ptr,
                                                        control_points_ptr,
                                                        weights_ptr);
    para_dim_ = c_spline_->SplinepyParaDim();
    dim_ = c_spline_->SplinepyDim();
  }

  /// will throw if c_spline_ is not initialized.
  /// use this for runtime core calls
  CoreSpline_& Core() {
    if (!c_spline_) {
      splinepy::utils::PrintAndThrowError(
          "Core spline does not exist.",
          "Please first intialize core spline.");
    }

    return c_spline_;
  }

  const CoreSpline_& Core() const {
    if (!c_spline_) {
      splinepy::utils::PrintAndThrowError(
          "Core spline does not exist.",
          "Please first intialize core spline.");
    }

    return c_spline_;
  }

  std::string WhatAmI() const { return Core()->SplinepyWhatAmI(); }
  std::string Name() const { return Core()->SplinepySplineName(); }
  bool HasKnotVectors() const { return Core()->SplinepyHasKnotVectors(); }
  bool IsRational() const { return Core()->SplinepyIsRational(); }

  /// Returns currunt properties of core spline
  /// similar to update_p
  py::dict CurrentCoreProperties() const {
    py::dict dict_spline;

    // prepare property arrays
    // first, degrees and control_points
    py::array_t<int> degrees(para_dim_);
    int* degrees_ptr = static_cast<int*>(degrees.request().ptr);
    const int ncps = Core()->SplinepyNumberOfControlPoints();
    py::array_t<double> control_points(ncps * dim_);
    double* control_points_ptr =
        static_cast<double*>(control_points.request().ptr);

    // and maybes
    std::vector<std::vector<double>> knot_vectors;
    std::vector<std::vector<double>>* knot_vectors_ptr = nullptr;
    py::array_t<double> weights;
    double* weights_ptr = nullptr;

    if (Core()->SplinepyHasKnotVectors()) {
      knot_vectors_ptr = &knot_vectors;
    }
    if (Core()->SplinepyIsRational()) {
      weights = py::array_t<double>(ncps);
      weights_ptr = static_cast<double*>(weights.request().ptr);
    }

    Core()->SplinepyCurrentProperties(degrees_ptr,
                                      knot_vectors_ptr,
                                      control_points_ptr,
                                      weights_ptr);

    // process
    dict_spline["degrees"] = degrees;
    control_points.resize({ncps, dim_});
    dict_spline["control_points"] = control_points;

    // process maybes
    if (Core()->SplinepyHasKnotVectors()) {
      py::list kvs;
      for (const auto& knot_vector : knot_vectors) {
        const std::size_t kvsize = knot_vector.size();
        py::array_t<double> kv(kvsize);
        double* kv_ptr = static_cast<double*>(kv.request().ptr);
        std::copy_n(knot_vector.begin(), kvsize, kv_ptr);
        kvs.append(kv);
      }
      dict_spline["knot_vectors"] = kvs;
    }
    if (Core()->SplinepyIsRational()) {
      weights.resize({ncps, 1});
      dict_spline["weights"] = weights;
    }

    return dict_spline;
  }

  /// AABB of spline parametric space
  py::array_t<double> ParametricBounds() const {
    // prepare output - [[lower_bound], [upper_bound]]
    py::array_t<double> pbounds(2 * para_dim_);
    double* pbounds_ptr = static_cast<double*>(pbounds.request().ptr);

    Core()->SplinepyParametricBounds(pbounds_ptr);

    pbounds.resize({2, para_dim_});
    return pbounds;
  }

  py::array_t<int> ControlMeshResolutions() const {
    // prepare output
    py::array_t<int> cmr(para_dim_);
    int* cmr_ptr = static_cast<int*>(cmr.request().ptr);

    Core()->SplinepyControlMeshResolutions(cmr_ptr);

    return cmr;
  }

  /// @brief Evaluate Splines at boundary face centers
  /// @return numpy array with results
  py::array_t<double> EvaluateBoundaryCenters() const {
    // prepare output
    const int n_faces = 2 * para_dim_;
    py::array_t<double> queries(n_faces * para_dim_);
    double* queries_ptr = static_cast<double*>(queries.request().ptr);
    py::array_t<double> face_centers(n_faces * dim_);
    double* face_centers_ptr = static_cast<double*>(face_centers.request().ptr);

    for (int i{}; i < para_dim_; i++) {
      for (int j{}; j < para_dim_; j++) {
        if (i == j) {
          queries_ptr[2 * i * para_dim_ + j] = 0.;
          queries_ptr[(2 * i + 1) * para_dim_ + j] = 1.;
        } else {
          queries_ptr[2 * i * para_dim_ + j] = .5;
          queries_ptr[(2 * i + 1) * para_dim_ + j] = .5;
        }
      }
    }
    for (int i{}; i < n_faces; ++i) {
      Core()->SplinepyEvaluate(&queries_ptr[i * para_dim_],
                               &face_centers_ptr[i * dim_]);
    }

    face_centers.resize({n_faces, dim_});
    return face_centers;
  }

  py::array_t<double> Evaluate(py::array_t<double> queries,
                               int nthreads) const {
    CheckPyArrayShape(queries, {-1, para_dim_}, true);

    // prepare output
    const int n_queries = queries.shape(0);
    py::array_t<double> evaluated(n_queries * dim_);
    double* evaluated_ptr = static_cast<double*>(evaluated.request().ptr);

    // prepare vectorized evaluate queries
    double* queries_ptr = static_cast<double*>(queries.request().ptr);
    auto evaluate = [&](int begin, int end) {
      for (int i{begin}; i < end; ++i) {
        Core()->SplinepyEvaluate(&queries_ptr[i * para_dim_],
                                 &evaluated_ptr[i * dim_]);
      }
    };

    splinepy::utils::NThreadExecution(evaluate, n_queries, nthreads);

    evaluated.resize({n_queries, dim_});
    return evaluated;
  }

  /// Sample wraps evaluate to allow nthread executions
  /// Requires SplinepyParametricBounds
  py::array_t<double> Sample(py::array_t<int> resolutions, int nthreads) const {
    CheckPyArraySize(resolutions, para_dim_, true);

    // get sampling bounds
    std::vector<double> bounds(para_dim_ * 2);
    double* bounds_ptr = bounds.data();
    Core()->SplinepyParametricBounds(bounds_ptr);
    // prepare sampler
    auto grid = splinepy::utils::CStyleArrayPointerGridPoints(
        para_dim_,
        bounds_ptr,
        static_cast<int*>(resolutions.request().ptr));
    // prepare output
    const int n_sampled = grid.Size();
    py::array_t<double> sampled(n_sampled * dim_);
    double* sampled_ptr = static_cast<double*>(sampled.request().ptr);

    // wrap evaluate
    auto sample = [&](int begin, int end) {
      std::vector<double> query(para_dim_);
      double* query_ptr = query.data();
      for (int i{begin}; i < end; ++i) {
        grid.IdToGridPoint(i, query_ptr);
        Core()->SplinepyEvaluate(&query[0], &sampled_ptr[i * dim_]);
      }
    };

    splinepy::utils::NThreadExecution(sample, n_sampled, nthreads);

    sampled.resize({n_sampled, dim_});
    return sampled;
  }

  /// spline derivatives
  py::array_t<double> Derivative(py::array_t<double> queries,
                                 py::array_t<int> orders,
                                 int nthreads) const {
    // process input
    CheckPyArrayShape(queries, {-1, para_dim_}, true);
    const int n_queries = queries.shape(0);
    int constant_orders_factor = 0;
    if (!CheckPyArraySize(orders, para_dim_, false)) {
      if (!CheckPyArrayShape(orders, {n_queries, para_dim_}, false)) {
        splinepy::utils::PrintAndThrowError(
            "Derivative-query-orders (orders) must either have same size as",
            "spline's parametric dimension or same shape as queries.",
            "Expected: size{",
            para_dim_,
            "} or shape{",
            n_queries,
            ",",
            para_dim_,
            "}",
            "Given: size{",
            orders.size(),
            "}, approx. shape {(size / para_dim), para_dim} = {",
            orders.size() / para_dim_,
            ",",
            para_dim_,
            "}");
      } else {
        // orders shape matches query shape.
        constant_orders_factor = 1;
      }
    }

    // prepare output
    py::array_t<double> derived(n_queries * dim_);
    double* derived_ptr = static_cast<double*>(derived.request().ptr);

    // prepare lambda for nthread exe
    double* queries_ptr = static_cast<double*>(queries.request().ptr);
    int* orders_ptr = static_cast<int*>(orders.request().ptr);
    auto derive = [&](int begin, int end) {
      for (int i{begin}; i < end; ++i) {
        Core()->SplinepyDerivative(
            &queries_ptr[i * para_dim_],
            &orders_ptr[constant_orders_factor * i * para_dim_],
            &derived_ptr[i * dim_]);
      }
    };

    splinepy::utils::NThreadExecution(derive, n_queries, nthreads);

    derived.resize({n_queries, dim_});
    return derived;
  }

  /// Basis function values and support id
  py::tuple BasisAndSupport(py::array_t<double> queries, int nthreads) {
    CheckPyArrayShape(queries, {-1, para_dim_}, true);
    const int n_queries = queries.shape(0);

    // prepare results
    const int n_support = Core()->SplinepyNumberOfSupports();
    py::array_t<double> basis(n_queries * n_support);
    py::array_t<int> support(n_queries * n_support);

    // prepare_lambda for nthread exe
    double* queries_ptr = static_cast<double*>(queries.request().ptr);
    double* basis_ptr = static_cast<double*>(basis.request().ptr);
    int* support_ptr = static_cast<int*>(support.request().ptr);
    auto basis_support = [&](int begin, int end) {
      for (int i{begin}; i < end; ++i) {
        Core()->SplinepyBasisAndSupport(&queries_ptr[i * para_dim_],
                                        &basis_ptr[i * n_support],
                                        &support_ptr[i * n_support]);
      }
    };

    splinepy::utils::NThreadExecution(basis_support, n_queries, nthreads);

    basis.resize({n_queries, n_support});
    support.resize({n_queries, n_support});

    return py::make_tuple(basis, support);
  }

  /// Basis function values and support id
  py::tuple BasisDerivativeAndSupport(py::array_t<double> queries,
                                      py::array_t<int> orders,
                                      int nthreads) {
    CheckPyArrayShape(queries, {-1, para_dim_}, true);
    const int n_queries = queries.shape(0);

    int constant_orders_factor = 0;
    if (!CheckPyArraySize(orders, para_dim_, false)) {
      if (!CheckPyArrayShape(orders, {n_queries, para_dim_}, false)) {
        splinepy::utils::PrintAndThrowError(
            "Derivative-query-orders (orders) must either have same size as",
            "spline's parametric dimension or same shape as queries.",
            "Expected: size{",
            para_dim_,
            "} or shape{",
            n_queries,
            ",",
            para_dim_,
            "}",
            "Given: size{",
            orders.size(),
            "}, approx. shape {(size / para_dim), para_dim} = {",
            orders.size() / para_dim_,
            ",",
            para_dim_,
            "}");
      } else {
        constant_orders_factor = 1;
      }
    }

    // prepare results
    const int n_support = Core()->SplinepyNumberOfSupports();
    py::array_t<double> basis_der(n_queries * n_support);
    py::array_t<int> support(n_queries * n_support);

    // prepare_lambda for nthread exe
    double* queries_ptr = static_cast<double*>(queries.request().ptr);
    int* orders_ptr = static_cast<int*>(orders.request().ptr);
    double* basis_der_ptr = static_cast<double*>(basis_der.request().ptr);
    int* support_ptr = static_cast<int*>(support.request().ptr);
    auto basis_der_support = [&](int begin, int end) {
      for (int i{begin}; i < end; ++i) {
        Core()->SplinepyBasisDerivativeAndSupport(
            &queries_ptr[i * para_dim_],
            &orders_ptr[constant_orders_factor * i * para_dim_],
            &basis_der_ptr[i * n_support],
            &support_ptr[i * n_support]);
      }
    };

    splinepy::utils::NThreadExecution(basis_der_support, n_queries, nthreads);

    basis_der.resize({n_queries, n_support});
    support.resize({n_queries, n_support});

    return py::make_tuple(basis_der, support);
  }

  /// Proximity query (verbose)
  py::tuple Proximities(py::array_t<double> queries,
                        py::array_t<int> initial_guess_sample_resolutions,
                        double tolerance,
                        int max_iterations,
                        bool aggresive_search_bounds,
                        int nthreads) {
    CheckPyArrayShape(queries, {-1, dim_}, true);
    CheckPyArraySize(initial_guess_sample_resolutions, para_dim_);

    const int n_queries = queries.shape(0);
    const int pd = para_dim_ * dim_;
    const int ppd = para_dim_ * pd;

    // prepare results
    py::array_t<double> para_coord(n_queries * para_dim_);
    py::array_t<double> phys_coord(n_queries * dim_);
    py::array_t<double> phys_diff(n_queries * dim_);
    py::array_t<double> distance(n_queries);
    py::array_t<double> convergence_norm(n_queries);
    py::array_t<double> first_derivatives(n_queries * pd);
    py::array_t<double> second_derivatives(n_queries * ppd);

    // prepare lambda for nthread exe
    double* queries_ptr = static_cast<double*>(queries.request().ptr);
    double* para_coord_ptr = static_cast<double*>(para_coord.request().ptr);
    double* phys_coord_ptr = static_cast<double*>(phys_coord.request().ptr);
    double* phys_diff_ptr = static_cast<double*>(phys_diff.request().ptr);
    double* distance_ptr = static_cast<double*>(distance.request().ptr);
    double* convergence_norm_ptr =
        static_cast<double*>(convergence_norm.request().ptr);
    double* first_derivatives_ptr =
        static_cast<double*>(first_derivatives.request().ptr);
    double* second_derivatives_ptr =
        static_cast<double*>(second_derivatives.request().ptr);
    auto proximities = [&](int begin, int end) {
      for (int i{begin}; i < end; ++i) {
        Core()->SplinepyVerboseProximity(&queries_ptr[i * dim_],
                                         tolerance,
                                         max_iterations,
                                         aggresive_search_bounds,
                                         &para_coord_ptr[i * para_dim_],
                                         &phys_coord_ptr[i * dim_],
                                         &phys_diff_ptr[i * dim_],
                                         distance_ptr[i],
                                         convergence_norm_ptr[i],
                                         &first_derivatives_ptr[i * pd],
                                         &second_derivatives_ptr[i * ppd]);
      }
    };

    // make sure kdtree is planted before query.
    //
    // there are two cases, where tree will not be built:
    // 1. same resolution as last proximity call - core spline will check
    // 2. any negative entry in resolutions - checked here
    //
    // there is one cases, where runtime_error will be thrown:
    // 1. entry smaller than 2 - checked here
    //
    // allow us to use abbreviation here.
    int* igsr_ptr =
        static_cast<int*>(initial_guess_sample_resolutions.request().ptr);
    bool plant_kdtree = true;
    for (int i{}; i < para_dim_; ++i) {
      const int& res = igsr_ptr[i];
      if (res < 0) {
        plant_kdtree = false;
        break;
      }
      if (res < 2) {
        splinepy::utils::PrintAndThrowError(
            "positive entries for initial_guess_sample_resolutions",
            "shouldn't be smaller than 2.",
            "Entry number: [",
            i,
            "] is",
            res);
      }
    }
    // yes, we could've built an input and called the function directly,
    // but, we will stick with calling interface functions
    if (plant_kdtree) {
      Core()->SplinepyPlantNewKdTreeForProximity(igsr_ptr, nthreads);
    }

    splinepy::utils::NThreadExecution(proximities, n_queries, nthreads);

    para_coord.resize({n_queries, para_dim_});
    phys_coord.resize({n_queries, para_dim_});
    phys_diff.resize({n_queries, dim_});
    distance.resize({n_queries, 1});
    convergence_norm.resize({n_queries, 1});
    first_derivatives.resize({n_queries, para_dim_, dim_});
    second_derivatives.resize({n_queries, para_dim_, para_dim_, dim_});

    return py::make_tuple(para_coord,
                          phys_coord,
                          phys_diff,
                          distance,
                          convergence_norm,
                          first_derivatives,
                          second_derivatives);
  }

  /// (multiple) Degree elevation
  void ElevateDegrees(py::array_t<int> para_dims) {
    int* para_dims_ptr = static_cast<int*>(para_dims.request().ptr);
    const int n_request = para_dims.size();
    for (int i{}; i < n_request; ++i) {
      const int& p_dim = para_dims_ptr[i];
      if (!(p_dim < para_dim_) || p_dim < 0) {
        splinepy::utils::PrintAndThrowError(
            p_dim,
            "is invalid parametric dimension for degree elevation.");
      }
      Core()->SplinepyElevateDegree(p_dim);
    }
  }

  /// (multiple) Degree Reduction
  /// returns a list of reduction result (bool)
  py::list ReduceDegrees(py::array_t<int> para_dims, double tolerance) {
    int* para_dims_ptr = static_cast<int*>(para_dims.request().ptr);
    const int n_request = para_dims.size();

    py::list successful;
    for (int i{}; i < n_request; ++i) {
      const int& p_dim = para_dims_ptr[i];
      if (!(p_dim < para_dim_) || p_dim < 0) {
        splinepy::utils::PrintAndThrowError(
            p_dim,
            "is invalid parametric dimension for degree reduction.");
      }
      successful.append(Core()->SplinepyReduceDegree(p_dim, tolerance));
    }

    return successful;
  }
};

/* operations for certain splines */

/// (multiple) knot insertion, single dimension
inline void
InsertKnots(PySpline& spline, int para_dim, py::array_t<double> knots) {
  double* knots_ptr = static_cast<double*>(knots.request().ptr);
  const int n_request = knots.size();

  for (int i{}; i < n_request; ++i) {
    spline.Core()->SplinepyInsertKnot(para_dim, knots_ptr[i]);
  }
}

/// (multiple) knot removal, single dimension
inline py::list RemoveKnots(PySpline& spline,
                            int para_dim,
                            py::array_t<double> knots,
                            double tolerance) {
  double* knots_ptr = static_cast<double*>(knots.request().ptr);
  const int n_request = knots.size();

  py::list successful;
  for (int i{}; i < n_request; ++i) {
    successful.append(
        spline.Core()->SplinepyRemoveKnot(para_dim, knots_ptr[i], tolerance));
  }

  return successful;
}

/// spline multiplication - currently only for bezier
inline PySpline Multiply(const PySpline& a, const PySpline& b) {
  // performs runtime checks and throws error
  return PySpline(a.Core()->SplinepyMultiply(b.Core()));
}

/// spline addition - currently only for bezier
inline PySpline Add(const PySpline& a, const PySpline& b) {
  // performs runtime checks and throws error
  return PySpline(a.Core()->SplinepyAdd(b.Core()));
}

/// spline composition - currently only for bezier
inline PySpline Compose(const PySpline& inner, const PySpline& outer) {
  // performs runtime checks and throws error
  return PySpline(inner.Core()->SplinepyCompose(outer.Core()));
}

/// spline derivative spline
inline PySpline DerivativeSpline(const PySpline& spline,
                                 py::array_t<int> orders) {
  CheckPyArraySize(orders, spline.para_dim_);

  int* orders_ptr = static_cast<int*>(orders.request().ptr);
  return PySpline(spline.Core()->SplinepyDerivativeSpline(orders_ptr));
}

/// spline split - returns py::list of PySplines
inline py::list
Split(const PySpline& spline, int p_dim, py::array_t<double> locations) {
  // make sure they are sorted
  std::vector<double> locs(locations.size());
  std::copy_n(static_cast<double*>(locations.request().ptr),
              locations.size(),
              locs.begin());
  // sort
  std::sort(locs.begin(), locs.end());

  // split and append
  py::list splitted;
  // very first
  auto tmp_splitted = spline.Core()->SplinepySplit(p_dim, locs[0]);
  splitted.append(PySpline(tmp_splitted[0]));
  for (std::size_t i{1}; i < locs.size(); ++i) {
    const double locs_with_offset =
        (locs[i] - locs[i - 1]) / (1. - locs[i - 1]);
    tmp_splitted = tmp_splitted[1]->SplinepySplit(p_dim, locs_with_offset);
    splitted.append(PySpline(tmp_splitted[0]));
  }
  // very last
  splitted.append(PySpline(tmp_splitted[1]));

  return splitted;
}

/// bezier patch extraction.
inline py::list ExtractBezierPatches(const PySpline& spline) {
  const auto sp_patches = spline.Core()->SplinepyExtractBezierPatches();
  py::list patches;
  for (const auto& p : sp_patches) {
    patches.append(PySpline(p));
  }
  return patches;
}

/// boundary spline extraction
inline PySpline ExtractBoundaries(const PySpline& spline,
                                  const py::array_t<int>& boundary_ids) {
  // Init return value
  py::list boundary_splines{};
  const int n_boundaries = boundary_ids.size();
  int* bid_ptr = static_cast<int*>(boundary_ids.request().ptr);
  if (boundary_ids.size() == 0) {
    for (int i{}; i < para_dim_ * 2; ++i) {
      boudnary_splines.append(
          PySpline(spline.Core()->SplinepyExtractBoundary(i)));
    }
  } else {
    const int max_bid = para_dim_ * 2 - 1;
    for (int i{}; i < n_boundaries; ++i) {
      const int& bid = bid_ptr[i];
      if (bid < 0 || bid > max_bid) {
        splinepy::utils::PrintAndThrowError("Requested Boundary ID :",
                                            bid,
                                            "exceeds admissible range.");
      }
      boudnary_splines.append(
          PySpline(spline.Core()->SplinepyExtractBoundary(bid_ptr[i])));
    }
  }

  return boundary_splines;
}

/// boundary spline extraction from axis and extreme
inline PySpline ExtractBoundaryFromAxisAndExtrema(const PySpline& spline,
                                                  const int& axis,
                                                  const int& extreme) {
  // Determine corresponding ID
  const int boundary_id = (extreme > 0) ? 2 * axis + 1 : 2 * axis;

  // Extract boundary
  return PySpline(spline.Core()->SplinepyExtractBoundary(boundary_id));
}

/// extract a single physical dimension from a spline
inline PySpline ExtractDim(const PySpline& spline, int phys_dim) {
  return PySpline(spline.Core()->SplinepyExtractDim(phys_dim));
}

/// composition derivative
inline PySpline CompositionDerivative(const PySpline& outer,
                                      const PySpline& inner,
                                      const PySpline& inner_derivative) {
  return PySpline(
      outer.Core()->SplinepyCompositionDerivative(inner.Core(),
                                                  inner_derivative.Core()));
}

/// returns a spline with knot vectors.
/// if the spline already has knots, it returns the same spline
/// else, returns a same spline with knots
inline PySpline SameSplineWithKnotVectors(PySpline& spline) {
  // early exit if the spline has knot vectors already
  if (spline.HasKnotVectors()) {
    return spline;
  }

  py::dict props = spline.CurrentCoreProperties();

  // based on degrees, generate knot vectors
  py::array_t<int> degrees = py::cast<py::array_t<int>>(props["degrees"]);
  int* degrees_ptr = static_cast<int*>(degrees.request().ptr);
  py::list kvs;
  for (int i{}; i < degrees.size(); ++i) {
    // prepare kv and get required number of repeating knots
    const int n_knot_repeat = degrees_ptr[i] + 1;
    py::array_t<double> kv(n_knot_repeat * 2);
    double* kv_ptr = static_cast<double*>(kv.request().ptr);

    // defined knots with 0 and 1
    for (int j{}; j < n_knot_repeat; ++j) {
      kv_ptr[j] = 0.;
      kv_ptr[n_knot_repeat + j] = 1.;
    }

    kvs.append(kv);
  }

  // update knot vectors to dict spline
  props["knot_vectors"] = kvs;

  return PySpline(props);
}

/// returns core spline's ptr address
inline intptr_t CoreId(const PySpline& spline) {
  return reinterpret_cast<intptr_t>(spline.Core().get());
}

/// reference count of core spline
inline int CoreRefCount(const PySpline& spline) {
  return spline.Core().use_count();
}

/// have core? A non error raising checker
inline bool HaveCore(const PySpline& spline) {
  return (spline.c_spline_) ? true : false;
}

/// Overwrite core with a nullptr and assign neg values to dims
inline void AnnulCore(PySpline& spline) {
  spline.c_spline_ = nullptr;
  spline.para_dim_ = -1;
  spline.dim_ = -1;
}

/// Internal use only
/// Extract CoreSpline_s from list of PySplines
inline std::vector<PySpline::CoreSpline_>
ListOfPySplinesToVectorOfCoreSplines(py::list pysplines) {
  // prepare return obj
  std::vector<PySpline::CoreSpline_> core_splines;
  core_splines.reserve(pysplines.size());

  // loop and append
  for (py::handle pys : pysplines) {
    core_splines.emplace_back(py::cast<PySpline>(pys).Core());
  }

  return core_splines;
}

inline void add_spline_pyclass(py::module& m, const char* class_name) {
  py::class_<splinepy::py::PySpline> klasse(m, class_name);

  klasse.def(py::init<>())
      .def(py::init<py::kwargs>()) // doc here?
      .def(py::init<splinepy::py::PySpline&>())
      .def("new_core", &splinepy::py::PySpline::NewCore)
      .def_readwrite("_data", &splinepy::py::PySpline::data_)
      .def_readonly("para_dim", &splinepy::py::PySpline::para_dim_)
      .def_readonly("dim", &splinepy::py::PySpline::dim_)
      .def_property_readonly("whatami", &splinepy::py::PySpline::WhatAmI)
      .def_property_readonly("name", &splinepy::py::PySpline::Name)
      .def_property_readonly("has_knot_vectors",
                             &splinepy::py::PySpline::HasKnotVectors)
      .def_property_readonly("is_rational", &splinepy::py::PySpline::IsRational)
      .def_property_readonly("parametric_bounds",
                             &splinepy::py::PySpline::ParametricBounds)
      .def_property_readonly("control_mesh_resolutions",
                             &splinepy::py::PySpline::ControlMeshResolutions)
      .def_property_readonly("evaluate_boundary_centers",
                             &splinepy::py::PySpline::EvaluateBoundaryCenters)
      .def("current_core_properties",
           &splinepy::py::PySpline::CurrentCoreProperties)
      .def("evaluate",
           &splinepy::py::PySpline::Evaluate,
           py::arg("queries"),
           py::arg("nthreads") = 1)
      .def("sample",
           &splinepy::py::PySpline::Sample,
           py::arg("resolutions"),
           py::arg("nthreads") = 1)
      .def("derivative",
           &splinepy::py::PySpline::Derivative,
           py::arg("queries"),
           py::arg("orders"),
           py::arg("nthreads") = 1)
      .def("basis_and_support",
           &splinepy::py::PySpline::BasisAndSupport,
           py::arg("queries"),
           py::arg("nthreads") = 1)
      .def("basis_derivative_and_support",
           &splinepy::py::PySpline::BasisDerivativeAndSupport,
           py::arg("queries"),
           py::arg("orders"),
           py::arg("nthreads") = 1)
      .def("proximities",
           &splinepy::py::PySpline::Proximities,
           py::arg("queries"),
           py::arg("initial_guess_sample_resolutions"),
           py::arg("tolerance"),
           py::arg("max_iterations") = -1,
           py::arg("aggressive_search_bounds") = false,
           py::arg("nthreads") = 1)
      .def("elevate_degrees",
           &splinepy::py::PySpline::ElevateDegrees,
           py::arg("para_dims"))
      .def("reduce_degrees",
           &splinepy::py::PySpline::ReduceDegrees,
           py::arg("para_dims"),
           py::arg("tolerance"))
      .def(py::pickle(
          [](splinepy::py::PySpline& spl) {
            return py::make_tuple(spl.CurrentCoreProperties(), spl.data_);
          },
          [](py::tuple t) {
            if (t.size() != 2) {
              splinepy::utils::PrintAndThrowError("Invalid PySpline state.");
            }

            py::kwargs properties = t[0].cast<py::kwargs>(); // init
            PySpline spl(properties);
            spl.data_ = t[1].cast<py::dict>(); // saved data

            return spl;
          }));

  m.def("insert_knots",
        &splinepy::py::InsertKnots,
        py::arg("spline"),
        py::arg("para_dim"),
        py::arg("knots"));
  m.def("remove_knots",
        &splinepy::py::RemoveKnots,
        py::arg("spline"),
        py::arg("para_dim"),
        py::arg("knots"),
        py::arg("tolerance"));
  m.def("multiply", &splinepy::py::Multiply, py::arg("a"), py::arg("b"));
  m.def("add", &splinepy::py::Add, py::arg("a"), py::arg("b"));
  m.def("compose", &splinepy::py::Compose, py::arg("outer"), py::arg("inner"));
  m.def("derivative_spline",
        &splinepy::py::DerivativeSpline,
        py::arg("spline"),
        py::arg("orders"));
  m.def("split",
        &splinepy::py::Split,
        py::arg("spline"),
        py::arg("para_dim"),
        py::arg("locations"));
  m.def("extract_bezier_patches",
        &splinepy::py::ExtractBezierPatches,
        py::arg("spline"));
  m.def("extract_boundary",
        &splinepy::py::ExtractBoundary,
        py::arg("spline"),
        py::arg("boundary_id"));
  m.def("extract_dim",
        &splinepy::py::ExtractDim,
        py::arg("spline"),
        py::arg("phys_dim"));
  m.def("composition_derivative",
        &splinepy::py::CompositionDerivative,
        py::arg("outer"),
        py::arg("inner"),
        py::arg("inner_derivative"));
  m.def("same_spline_with_knot_vectors",
        &splinepy::py::SameSplineWithKnotVectors,
        py::arg("spline"));
  m.def("core_id", &splinepy::py::CoreId, py::arg("spline"));
  m.def("core_ref_count", &splinepy::py::CoreRefCount, py::arg("spline"));
  m.def("have_core", &splinepy::py::HaveCore, py::arg("spline"));
  m.def("annul_core", &splinepy::py::AnnulCore, py::arg("spline"));
  ;
}

} // namespace splinepy::py
