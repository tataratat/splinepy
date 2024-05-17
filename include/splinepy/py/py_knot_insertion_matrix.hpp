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
 * @param new_knots new knots
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
