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

template<typename ValueType>
static bool CheckPyArrayShape(const py::array_t<ValueType> arr,
                              const std::vector<int>& shape,
                              const bool throw_ = true) {
  const std::size_t expected_dim = shape.size();
  if (expected_dim != arr.ndim()) {
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

  CoreSpline_ c_spline_;
  int para_dim_;
  int dim_;
  bool is_rational_;
  bool has_knot_vectors_;

  // ctor
  PySpline() = default;
  PySpline(const py::kwargs& kwargs) { NewCore(kwargs); }
  PySpline(const CoreSpline_& another_core) : c_spline_(another_core) {
    para_dim_ = c_spline_->SplinepyParaDim();
    dim_ = c_spline_->SplinepyDim();
    is_rational_ = c_spline_->SplinepyIsRational();
    has_knot_vectors_ = c_spline_->SplinepyHasKnotVectors();
  }

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
    has_knot_vectors_ = false;

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
      has_knot_vectors_ = true;
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

    if (has_knot_vectors_) {
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
    if (has_knot_vectors_) {
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
    CheckPyArrayShape(queries, {-1, para_dim_}, true);

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
    CheckPyArraySize(resolutions, para_dim_, true);

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
    CheckPyArrayShape(queries, {-1, para_dim_}, true);
    const int n_queries = queries.shape(0);
    int constant_orders_factor = 0;
    if (!CheckPyArraySize(orders, para_dim_, false)) {
      if (!CheckPyArrayShape(orders, {n_queries, para_dim_}, false)) {
        constant_orders_factor = 1;
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

  /// Basis function values and support id
  py::tuple BasisAndSupport(py::array_t<double> queries, int nthreads) {
    CheckPyArrayShape(queries, {-1, para_dim_}, true);
    const int n_queries = queries.shape(0);

    // prepare results
    const int n_support = c_spline_->SplinepyNumberOfSupports();
    py::array_t<double> basis(n_queries * n_support);
    py::array_t<int> support(n_queries * n_support);

    // prepare_lambda for nthread exe
    double* queries_ptr = static_cast<double*>(queries.request().ptr);
    double* basis_ptr = static_cast<double*>(basis.request().ptr);
    int* support_ptr = static_cast<int*>(support.request().ptr);
    auto basis_support = [&](int begin, int end) {
      for (int i{begin}; i < end; ++i) {
        c_spline_->SplinepyBasisAndSupport(&queries_ptr[i * para_dim_],
                                           &basis_ptr[i * n_support],
                                           &support_ptr[i * n_support]);
      }
    };

    splinepy::utils::NThreadExecution(basis_support, n_queries, nthreads);

    basis.resize({n_queries, n_support});
    support.resize({n_queries, n_support});

    return py::make_tuple(basis, support);
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
        c_spline_->SplinepyVerboseProximity(&queries_ptr[i * dim_],
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
      c_spline_->SplinepyPlantNewKdTreeForProximity(igsr_ptr, nthreads);
    }

    splinepy::utils::NThreadExecution(proximities, n_queries, nthreads);

    para_coord.resize({n_queries, para_dim_});
    phys_coord.resize({n_queries, para_dim_});
    phys_diff.resize({n_queries, para_dim_});
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

/// spline multiplication - currently only for bezier
PySpline Multiply(const PySpline& a, const PySpline& b) {
  // performs runtime checks and throws error
  return PySpline(a.c_spline_->SplinepyMultiply(b.c_spline_));
}

/// spline addition - currently only for bezier
PySpline Add(const PySpline& a, const PySpline& b) {
  // performs runtime checks and throws error
  return PySpline(a.c_spline_->SplinepyAdd(b.c_spline_));
}

/// spline composition - currently only for bezier
PySpline Compose(const PySpline& inner, const PySpline& outer) {
  // performs runtime checks and throws error
  return PySpline(inner.c_spline_->SplinepyCompose(outer.c_spline_));
}

/// spline derivative spline
PySpline DerivativeSpline(const PySpline& spline, py::array_t<int> orders) {
  CheckPyArraySize(orders, spline.para_dim_);

  int* orders_ptr = static_cast<int*>(orders.request().ptr);
  return PySpline(spline.c_spline_->SplinepyDerivativeSpline(orders_ptr));
}

/// spline split - returns py::list of PySplines
py::list
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
  auto tmp_splitted = spline.c_spline_->SplinepySplit(p_dim, locs[0]);
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
py::list ExtractBezierPatches(const PySpline& spline) {
  const auto sp_patches = spline.c_spline_->SplinepyExtractBezierPatches();
  py::list patches;
  for (const auto& p : sp_patches) {
    patches.append(PySpline(p));
  }
  return patches;
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
      .def("basis_and_support",
           &splinepy::py::PySpline::BasisAndSupport,
           py::arg("queries"),
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
  ;
}
