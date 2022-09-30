#include <memory>
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

/// True interface to python
class PySpline {
public:
  using CoreSpline_ = typename std::shared_ptr<splinepy::splines::SplinepyBase>;

  CoreSpline_ c_spline_;
  int para_dim_;
  int dim_;
  bool is_rational_;
  bool is_bspline_;

  // ctor
  PySpline() = default;
  PySpline(const py::kwargs& kwargs) { NewCore(kwargs); }

  /// Creates a corresponding spline based on kwargs
  /// similar to previous update_c()
  void NewCore(const py::kwargs& kwargs) {
    // parse kwargs
    double* degrees_ptr = nullptr;
    std::vector<std::vector<double>> knot_vectors;
    std::vector<std::vector<double>>* knot_vectors_ptr = nullptr;
    double* control_points_ptr = nullptr;
    double* weights_ptr = nullptr;
    int para_dim, dim;
    is_rational_ = false;
    is_bspline_ = false;

    // get degrees and set para_dim
    auto d_array = py::cast<py::array_t<double>>(kwargs["degrees"]);
    degrees_ptr = static_cast<double*>(d_array.request().ptr);
    para_dim = d_array.shape(0);
    knot_vectors.reserve(para_dim);

    // get cps and set dim
    auto cp_array = py::cast<py::array_t<double>>(kwargs["control_points"]);
    control_points_ptr = static_cast<double*>(cp_array.request().ptr);
    dim = cp_array.shape(1);

    // maybe, get knot_vectors
    if (kwargs.contains("knot_vectors")) {
      is_bspline_ = true;
      for (py::handle kv : kwargs["knot_vectors"]) {
        std::vector<double> knot_vector;
        // cast to array_t, as python side will try to use tracked array anyways
        auto kv_array = py::cast<py::array_t<double>>(kv);
        knot_vector.reserve(kv_array.size());
        for (py::handle k : kv) {
          knot_vector.push_back(k.cast<double>());
        }
        knot_vectors.push_back(std::move(knot_vector));
      }
      knot_vectors_ptr = &knot_vectors;
    }

    // maybe, get weights
    if (kwargs.contains("weights")) {
      is_rational_ = true;
      auto w_array = py::cast<py::array_t<double>>(kwargs["weights"]);
      weights_ptr = static_cast<double*>(w_array.request().ptr);
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

  std::string WhatAmI() const { return c_spline_->SplinepyWhatAmI(); }

  // Returns currunt properties of core spline
  // similar to update_p
  py::dict CurrentProperties() const {
    py::dict dict_spline;

    // prepare property arrays
    // first, degrees and control_points
    py::array_t<double> degrees(para_dim_);
    double* degrees_ptr = static_cast<double*>(degrees.request().ptr);
    const int ncps = c_spline_->SplinepyNumberOfControlPoints();
    py::array_t<double> control_points(ncps * dim_);
    double* control_points_ptr =
        static_cast<double*>(control_points.request().ptr);

    // and maybes
    std::vector<std::vector<double>> knot_vectors;
    std::vector<std::vector<double>>* knot_vectors_ptr = nullptr;
    py::array_t<double> weights;
    double* weights_ptr = nullptr;

    if (is_bspline_) {
      knot_vectors_ptr = &knot_vectors;
    }
    if (is_rational_) {
      weights = py::array_t<double>(ncps);
      weights_ptr = static_cast<double*>(weights.request().ptr);
    }

    c_spline_->SplinepyCurrentProperties(degrees_ptr,
                                         knot_vectors_ptr,
                                         control_points_ptr,
                                         weights_ptr);

    // process
    dict_spline["degrees"] = degrees;
    control_points.resize({ncps, dim_});
    dict_spline["control_points"] = control_points;

    // process maybes
    if (is_bspline_) {
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
    if (is_rational_) {
      weights.resize({ncps, 1});
      dict_spline["weights"] = weights;
    }

    return dict_spline;
  }

  py::array_t<double> Evaluate(py::array_t<double> queries,
                               int nthreads) const {
    // prepare output
    const int n_queries = queries.shape(0);
    py::array_t<double> evaluated(n_queries * dim_);
    double* evaluated_ptr = static_cast<double*>(evaluated.request().ptr);

    // prepare vectorized evaluate queries
    double* queries_ptr = static_cast<double*>(queries.request().ptr);
    auto evaluate = [&](int begin, int end) {
      for (int i{begin}; i < end; ++i) {
        c_spline_->SplinepyEvaluate(&queries_ptr[i * para_dim_],
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
    // get sampling bounds
    std::vector<double> bounds(para_dim_ * 2);
    double* bounds_ptr = bounds.data();
    c_spline_->SplinepyParametricBounds(bounds_ptr);
    // prepare sampler
    auto grid = splinepy::utils::RawPtrGridPoints(
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
        c_spline_->SplinepyEvaluate(&query[0], &sampled_ptr[i * dim_]);
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
    const int n_queries = queries.shape(0);
    const int orders_ndim = orders.ndim();
    const int orders_len = orders.shape(0);
    int constant_orders_factor = 0;
    if (orders_ndim != 1) {
      if (orders_len == n_queries) {
        constant_orders_factor = 1;
      } else if (orders_len == 1) {
        // pass
      } else {
        splinepy::utils::PrintAndThrowError(
            "Length of derivative-query-orders (orders) must be either 1",
            "or same as the length of queries.",
            "Expected:",
            n_queries,
            "Given:",
            orders_len);
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
        c_spline_->SplinepyDerivative(
            &queries_ptr[i * para_dim_],
            &orders_ptr[constant_orders_factor * i * para_dim_],
            &derived_ptr[i * dim_]);
      }
    };

    splinepy::utils::NThreadExecution(derive, n_queries, nthreads);

    derived.resize({n_queries, dim_});
    return derived;
  }

  /// (multiple) Degree elevation
  void ElevateDegrees(py::array_t<int> para_dims) {
    int* para_dims_ptr = static_cast<int*>(para_dims.request().ptr);
    const int n_request = para_dims.size();
    for (int i{}; i < n_request; ++i) {
      c_spline_->SplinepyElevateDegree(para_dims_ptr[i]);
    }
  }

  /// (multiple) Degree Reduction
  /// returns a list of reduction result (bool)
  py::list ReduceDegrees(py::array_t<int> para_dims, double tolerance) {
    int* para_dims_ptr = static_cast<int*>(para_dims.request().ptr);
    const int n_request = para_dims.size();

    py::list successful;
    for (int i{}; i < n_request; ++i) {
      successful.append(
          c_spline_->SplinepyReduceDegree(para_dims_ptr[i], tolerance));
    }

    return successful;
  }
};

/* operations for certain splines */

/// (multiple) knot insertion, single dimension
void InsertKnots(PySpline& spline, int para_dim, py::array_t<double> knots) {
  double* knots_ptr = static_cast<double*>(knots.request().ptr);
  const int n_request = knots.size();

  for (int i{}; i < n_request; ++i) {
    spline.c_spline_->SplinepyInsertKnot(para_dim, knots_ptr[i]);
  }
}

/// (multiple) knot removal, single dimension
py::list RemoveKnots(PySpline& spline,
                     int para_dim,
                     py::array_t<double> knots,
                     double tolerance) {
  double* knots_ptr = static_cast<double*>(knots.request().ptr);
  const int n_request = knots.size();

  py::list successful;
  for (int i{}; i < n_request; ++i) {
    successful.append(spline.c_spline_->SplinepyRemoveKnot(para_dim,
                                                           knots_ptr[i],
                                                           tolerance));
  }

  return successful;
}

} // namespace splinepy::py

void add_spline_pyclass(py::module& m, const char* class_name) {
  py::class_<splinepy::py::PySpline> klasse(m, class_name);

  klasse.def(py::init<>())
      .def(py::init<py::kwargs>()) // doc here?
      .def_property_readonly("whatami", &splinepy::py::PySpline::WhatAmI)
      .def("current_properties", &splinepy::py::PySpline::CurrentProperties)
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
      .def("elevate_degrees",
           &splinepy::py::PySpline::ElevateDegrees,
           py::arg("para_dims"))
      .def("reduce_degrees",
           &splinepy::py::PySpline::ReduceDegrees,
           py::arg("para_dims"),
           py::arg("tolerance"));
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
        py::arg("tolernace"));
  ;
}
