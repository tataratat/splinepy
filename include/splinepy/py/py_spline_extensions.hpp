/// collection of functions that take spline as parameter
/// but either not applicable for all splines or only applicable to
/// specific use cases.
#pragma once

#include <memory>

// pybind11
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "splinepy/py/py_spline.hpp"
#include "splinepy/splines/null_spline.hpp"

namespace splinepy::py {

namespace py = pybind11;

/// (multiple) knot insertion, single dimension
py::list InsertKnots(std::shared_ptr<PySpline>& spline,
                     int para_dim,
                     py::array_t<double> knots);

/// (multiple) knot removal, single dimension
py::list RemoveKnots(std::shared_ptr<PySpline>& spline,
                     int para_dim,
                     py::array_t<double> knots,
                     double tolerance);

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

/// bezier patch extraction.
py::list ExtractBezierPatches(const std::shared_ptr<PySpline>& spline);

/// boundary spline extraction
py::list ExtractBoundaries(const std::shared_ptr<PySpline>& spline,
                           const py::array_t<int>& boundary_ids);

/// extract a single physical dimension from a spline
std::shared_ptr<PySpline> ExtractDim(const std::shared_ptr<PySpline>& spline,
                                     int phys_dim);

/// composition derivative
std::shared_ptr<PySpline>
CompositionDerivative(const std::shared_ptr<PySpline>& outer,
                      const std::shared_ptr<PySpline>& inner,
                      const std::shared_ptr<PySpline>& inner_derivative);

/// returns a spline with knot vectors.
/// if the spline already has knots, it returns the same spline
/// else, returns a same spline with knots
std::shared_ptr<PySpline>
SameSplineWithKnotVectors(std::shared_ptr<PySpline>& spline);

/// @brief Evaluate Splines at boundary face centers
/// @return numpy array with results
py::array_t<double> EvaluateBoundaryCenters(std::shared_ptr<PySpline>& spline);

/// returns core spline's ptr address
intptr_t CoreId(const std::shared_ptr<PySpline>& spline);

/// reference count of core spline
int CoreRefCount(const std::shared_ptr<PySpline>& spline);

/// have core? A non error raising checker
bool HasCore(const std::shared_ptr<PySpline>& spline);

/// Overwrite core with a nullptr and assign neg values to dims
void AnnulCore(std::shared_ptr<PySpline>& spline);

/// null spline creator
static inline std::shared_ptr<PySpline> CreateNullSpline(const int para_dim,
                                                         const int dim) {
  return std::make_shared<PySpline>(
      splinepy::splines::kNullSplineLookup[para_dim - 1][dim - 1]);
}

} // namespace splinepy::py
