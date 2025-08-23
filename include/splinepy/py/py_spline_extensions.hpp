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

/// returns a spline with knot vectors.
/// if the spline already has knots, it returns the same spline
/// else, returns a same spline with knots
std::shared_ptr<PySpline>
SameSplineWithKnotVectors(std::shared_ptr<PySpline>& spline);

/// @brief Evaluate Splines at boundary face centers
/// @return numpy array with results
py::array_t<double> EvaluateBoundaryCenters(std::shared_ptr<PySpline>& spline);

/// @brief Creates and returns NullSpline. Filler for multipatch
/// @param para_dim
/// @param dim
/// @return
static inline std::shared_ptr<PySpline> CreateNullSpline(const int para_dim,
                                                         const int dim) {
  return std::make_shared<PySpline>(
      splinepy::splines::kNullSplineLookup[para_dim - 1][dim - 1]);
}

} // namespace splinepy::py
