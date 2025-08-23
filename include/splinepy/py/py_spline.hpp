/*
MIT License

Copyright (c) 2021 Jaewook Lee

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#pragma once

// pybind11
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "splinepy/splines/splinepy_base.hpp"
#include "splinepy/splines/splinepy_bezier.hpp"
#include "splinepy/splines/splinepy_bspline.hpp"
#include "splinepy/splines/splinepy_rational.hpp"
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

class PySpline : public splinepy::splines::SplinepyBase {
public:
  using Base_ = splinepy::splines::SplinepyBase;

  /// AABB of spline parametric space
  virtual py::array_t<double> ParametricBounds() const;

  /// @brief Calculate Greville abscissae for Spline.
  virtual py::array_t<double> GrevilleAbscissae(const double) const;

  /// @brief Returns control mesh resolutions
  virtual py::array_t<int> ControlMeshResolutions() const;

  /// @brief Evaluate spline at query points
  /// @param queries Query points
  /// @param nthreads Number of threads to use
  virtual py::array_t<double> Evaluate(py::array_t<double> queries,
                                       int nthreads) const;

  /// Sample wraps evaluate to allow nthread executions
  /// Requires SplinepyParametricBounds
  virtual py::array_t<double> Sample(py::array_t<int> resolutions,
                                     int nthreads) const;

  /**
   * @brief Evaluate the Jacobian at certain positions
   *
   * @param queries position in the parametric space
   * @param nthreads number of threads for evaluation
   * @return py::array_t<double>
   */
  virtual py::array_t<double> Jacobian(const py::array_t<double> queries,
                                       const int nthreads) const;

  /// spline derivatives
  virtual py::array_t<double> Derivative(py::array_t<double> queries,
                                         py::array_t<int> orders,
                                         int nthreads) const;

  /// Basis support id
  virtual py::array_t<int> Support(py::array_t<double> queries,
                                   int nthreads) const;

  /// Basis function values
  virtual py::array_t<double> Basis(py::array_t<double> queries,
                                    int nthreads) const;

  /// Basis function values and support id
  virtual py::tuple BasisAndSupport(py::array_t<double> queries,
                                    int nthreads) const;

  /// @brief Get basis derivative
  /// @param queries Query points
  /// @param orders
  /// @param nthreads number of threads to use
  virtual py::array_t<double> BasisDerivative(py::array_t<double> queries,
                                              py::array_t<int> orders,
                                              int nthreads) const;

  /// Basis function values and support id
  virtual py::tuple BasisDerivativeAndSupport(py::array_t<double> queries,
                                              py::array_t<int> orders,
                                              int nthreads) const;

  /// Proximity query (verbose)
  virtual py::tuple
  Proximities(py::array_t<double> queries,
              py::array_t<int> initial_guess_sample_resolutions,
              const double tolerance,
              const int max_iterations,
              const bool aggresive_search_bounds,
              const int nthreads);

  /// (multiple) Degree elevation
  virtual void ElevateDegrees(py::array_t<int> para_dims);

  /// bezier patch extraction.
  py::list ExtractBezierPatches(const std::shared_ptr<PySpline>& spline);

  /// boundary spline extraction
  py::list ExtractBoundaries(const std::shared_ptr<PySpline>& spline,
                             const py::array_t<int>& boundary_ids);

  virtual py::array_t<int> Degrees() const;
  virtual py::object ControlPoints();

protected:
  /// @brief runtime registry for extension methods. Prio 0.
  static py::dict ext_methods_registry_;
  /// @brief runtime registry for helper classes. Prio 1.
  static py::dict ext_helper_registry_;
  /// @brief Subclassed (at python side) numpy class. Keeps track of inplace
  /// changes.
  py::object tracked_control_points_;
  /// @brief Subclassed (at python side) numpy class. Keeps track of inplace
  /// changes.
  py::object tracked_weights_;

  /// @brief Update backend control points
  virtual void SyncControlPoints();
  /// @brief Update tracked arrays (python)
  virtual void SyncTrackedArrays();
};
class PyRational : public splinepy::splines::SplinepyRational {
public:
  /// @brief Returns tracked weights.
  /// @return
  virtual py::object Weights();
};
class PyBezier : public PySpline, public splinepy::splines::SplinepyBezier {
public:
  static std::shared_ptr<PyBezier>
  Create(const py::array_t<int>& degrees,
         const py::array_t<double>& control_points);

  /// spline multiplication - currently only for bezier
  std::shared_ptr<PySpline> Multiply(const std::shared_ptr<PySpline>& a,
                                     const std::shared_ptr<PySpline>& b);

  /// spline addition - currently only for bezier
  std::shared_ptr<PySpline> Add(const std::shared_ptr<PySpline>& a,
                                const std::shared_ptr<PySpline>& b);

  /// spline composition - currently only for bezier
  std::shared_ptr<PySpline> Compose(const std::shared_ptr<PySpline>& outer,
                                    const std::shared_ptr<PySpline>& inner);

  /**
   * @brief Compute the sensitivities with respect to the control points of the
   * outer function given a composition
   *
   * @param inner Inner function (Bezier type)
   * @param outer Outer Function (Bezier type)
   * @return py::list list of Bezier splines representing the derivatives
   */
  py::list ComposeSensitivities(const std::shared_ptr<PySpline>& inner,
                                const std::shared_ptr<PySpline>& outer);

  /// spline derivative spline
  std::shared_ptr<PySpline>
  DerivativeSpline(const std::shared_ptr<PySpline>& spline,
                   py::array_t<int> orders);

  /// spline split - returns py::list of PySplines
  py::list Split(const std::shared_ptr<PySpline>& spline,
                 int p_dim,
                 py::array_t<double> locations);

  /// extract a single physical dimension from a spline
  std::shared_ptr<PySpline> ExtractDim(const std::shared_ptr<PySpline>& spline,
                                       int phys_dim);

  /// composition derivative
  std::shared_ptr<PySpline>
  CompositionDerivative(const std::shared_ptr<PySpline>& outer,
                        const std::shared_ptr<PySpline>& inner,
                        const std::shared_ptr<PySpline>& inner_derivative);
};

class PyRationalBezier : public PyBezier,
                         public PyRational,
                         public splinepy::splines::SplinepyBezier {
public:
  static std::shared_ptr<PyRationalBezier>
  Create(const py::array_t<int>& degrees,
         const py::array_t<double>& control_points,
         const py::array_t<double>& weights);

  virtual py::object Weights();
};

class PyBSpline : public PySpline, public splinepy::splines::SplinepyBSpline {
  static std::shared_ptr<PyBSpline>
  Create(const py::array_t<int>& degrees,
         const py::list& knot_vectors,
         const py::array_t<double>& control_points);

  bool HasKnotVectors() const;

  /// (multiple) Degree Reduction
  /// returns a list of reduction result (bool)
  py::list ReduceDegrees(py::array_t<int> para_dims, double tolerance);

  /// @brief Returns ParameterSpace. meant to be called library internally to
  /// prepare ParameterSpaceBase (forms knot_vectors)
  std::shared_ptr<bsplinelib::parameter_spaces::ParameterSpaceBase>
  ParameterSpace();

  /// (multiple) knot insertion, single dimension
  py::array_t<bool> InsertKnots(int para_dim, py::array_t<double> knots);

  /// (multiple) knot removal, single dimension
  py::list
  RemoveKnots(int para_dim, py::array_t<double> knots, double tolerance);
};
class PyNURBS : public PyBSpline,
                public PyRational public splinepy::splines::SplinepyBSpline {
  static std::shared_ptr<PyNURBS>
  Create(const py::array_t<int>& degrees,
         const py::list& knot_vectors,
         const py::array_t<double>& control_points,
         const py::array_t<double>& weights);

  virtual py::object Weights();
};
} // namespace splinepy::py
