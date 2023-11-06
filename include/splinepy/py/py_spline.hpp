#pragma once

// pybind11
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

// first four are required for Create* implementations
#include "splinepy/splines/splinepy_base.hpp"
#include "splinepy/utils/print.hpp"

namespace splinepy::py {

namespace py = pybind11;

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

  /// @brief Core spline
  CoreSpline_ c_spline_ = nullptr;

  // store very frequently used values - from python this will be readonly
  /// @brief Dimension of parameter space
  int para_dim_ = -1;
  /// @brief Dimension of physical space
  int dim_ = -1;

  /// @brief Place to store and access from both python and cpp side.
  py::dict data_;

  // ctor
  PySpline() = default;
  /// @brief Move constructor
  PySpline(PySpline&&) = default;
  /// @brief Constructor
  /// @param kwargs
  PySpline(const py::kwargs& kwargs) { NewCore(kwargs); }
  /// @brief Constructor
  /// @param another_core
  PySpline(const CoreSpline_& another_core) : c_spline_(another_core) {
    para_dim_ = c_spline_->SplinepyParaDim();
    dim_ = c_spline_->SplinepyDim();
  }
  /// @brief Constructor
  /// @param another_py_spline
  PySpline(PySpline& another_py_spline)
      : c_spline_(another_py_spline.Core()),
        para_dim_(another_py_spline.para_dim_),
        dim_(another_py_spline.dim_) {
    // nichts
  }
  /// @brief Constructor
  /// @param another_py_spline_ptr
  PySpline(const std::shared_ptr<PySpline>& another_py_spline_ptr)
      : PySpline(*another_py_spline_ptr) {}

  /// Creates a corresponding spline based on kwargs
  /// similar to previous update_c()
  /// Runs sanity checks on inputs
  void NewCore(const py::kwargs& kwargs);

  /// will throw if c_spline_ is not initialized.
  /// use this for runtime core calls
  CoreSpline_& Core();

  /// @brief Core spline
  const CoreSpline_& Core() const;

  /// @brief What am I?
  std::string WhatAmI() const { return Core()->SplinepyWhatAmI(); }
  /// @brief Get spline name
  std::string Name() const { return Core()->SplinepySplineName(); }
  /// @brief Returns true iff spline has knot vectors
  bool HasKnotVectors() const { return Core()->SplinepyHasKnotVectors(); }
  /// @brief Returns True iff spline is rational. NURBS is rational,
  /// for example.
  bool IsRational() const { return Core()->SplinepyIsRational(); }

  /// As knot vectors and control points / weights has a specific initialization
  /// routines, we provide a separate degree getter to avoid calling
  /// CurrentCoreProperties() for a full properties copy.
  py::array_t<int> CurrentCoreDegrees() const;

  /// Returns currunt properties of core spline
  /// similar to update_p
  py::dict CurrentCoreProperties() const;

  /// Returns coordinate pointers in a tuple. For rational splines,
  /// This will return control_point_pointers and weight_pointers
  /// For non-rational splines, only the former.
  py::tuple CoordinatePointers();

  /// Returns knot vector of given dimension. meant to be
  /// called library internally to prepare @property
  std::shared_ptr<bsplinelib::parameter_spaces::KnotVector>
  KnotVector(const int para_dim);

  /// AABB of spline parametric space
  py::array_t<double> ParametricBounds() const;

  /// @brief Calculate Greville abscissae for Spline
  py::list GrevilleAbscissae(const double) const;

  /// @brief Returns control mesh resolutions
  py::array_t<int> ControlMeshResolutions() const;

  /// @brief Evaluate spline at query points
  /// @param queries Query points
  /// @param nthreads Number of threads to use
  py::array_t<double> Evaluate(py::array_t<double> queries, int nthreads) const;

  /// Sample wraps evaluate to allow nthread executions
  /// Requires SplinepyParametricBounds
  py::array_t<double> Sample(py::array_t<int> resolutions, int nthreads) const;

  /**
   * @brief Evaluate the Jacobian at certain positions
   *
   * @param queries position in the parametric space
   * @param nthreads number of threads for evaluation
   * @return py::array_t<double>
   */
  py::array_t<double> Jacobian(const py::array_t<double> queries,
                               const int nthreads) const;

  /// spline derivatives
  py::array_t<double> Derivative(py::array_t<double> queries,
                                 py::array_t<int> orders,
                                 int nthreads) const;

  /// Basis support id
  py::array_t<int> Support(py::array_t<double> queries, int nthreads) const;

  /// Basis function values
  py::array_t<double> Basis(py::array_t<double> queries, int nthreads) const;

  /// Basis function values and support id
  py::tuple BasisAndSupport(py::array_t<double> queries, int nthreads) const;

  /// @brief Get basis derivative
  /// @param queries Query points
  /// @param orders
  /// @param nthreads number of threads to use
  py::array_t<double> BasisDerivative(py::array_t<double> queries,
                                      py::array_t<int> orders,
                                      int nthreads) const;

  /// Basis function values and support id
  py::tuple BasisDerivativeAndSupport(py::array_t<double> queries,
                                      py::array_t<int> orders,
                                      int nthreads) const;

  /// Proximity query (verbose)
  py::tuple Proximities(py::array_t<double> queries,
                        py::array_t<int> initial_guess_sample_resolutions,
                        double tolerance,
                        int max_iterations,
                        bool aggresive_search_bounds,
                        int nthreads);

  /// (multiple) Degree elevation
  void ElevateDegrees(py::array_t<int> para_dims);

  /// (multiple) Degree Reduction
  /// returns a list of reduction result (bool)
  py::list ReduceDegrees(py::array_t<int> para_dims, double tolerance);

  /// @brief returns current spline as package's derived spline types based on
  /// splinepy.settings.NAME_TO_TYPE
  /// @return
  py::object ToDerived();
};

} // namespace splinepy::py
