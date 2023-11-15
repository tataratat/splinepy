#pragma once

#include <algorithm>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
// for automatic conversion of py::list to std::vector
#include <pybind11/stl.h>

#include "splinepy/py/py_spline.hpp"

/// @brief
namespace splinepy::py {

namespace py = pybind11;

/**
 * @brief Create a local (dense) matrix that correlates two knot vectors
 *
 * @param old_kv Old knot vector
 * @param new_kv New knot vector
 * @param degree
 * @param tolerance tolerance (to identify multiple knots)
 * @return std::tuple<py::array_t<double>, std::vector<int>>
 */
std::tuple<py::array_t<double>, std::vector<int>>
ComputeKnotInsertionMatrixAndKnotSpan(const py::array_t<double>& old_kv,
                                      const py::array_t<double>& new_kv,
                                      const int degree,
                                      const double& tolerance);

/**
 * @brief  Compute the Global knot insertion matrix for all parametric
 * dimensions at once
 *
 * Currently we return a list of length (para_dim) containing a tuple of numpy
 * arrays, that can be used to instantiate scipy sparse matrices. This helps
 * avoid binding eigen as a separate library.
 *
 * @todo: replace calls with pybind/eigen
 *
 * @param old_kvs list of arrays, representing the individual knot vectors of
 * the start spline
 * @param degrees degrees along all parametric dimensions
 * @param parametric_dimension where knots are to be inserted
 * @param new_kvs new knots
 * @param tolerance tolerance for identifying individual knots
 * @return py::tuple to create new matrix
 */
py::tuple ComputeGlobalKnotInsertionMatrix(
    const std::vector<py::array_t<double>>& old_kvs,
    const py::array_t<int>& degrees,
    const int parametric_dimension,
    const py::array_t<double>& new_knots,
    const double& tolerance);

/**
 * @brief  Helper function to provide all necessary information to assemble and
 * compute the bezier extraction matrix, that links all control points to the
 * set of individual bezier patches
 *
 */
py::tuple
BezierExtractionMatrices(const std::shared_ptr<splinepy::py::PySpline>& spline,
                         const double& tolerance);
} // namespace splinepy::py
