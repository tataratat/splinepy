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

#include <pybind11/pybind11.h>

// core_spline
namespace splinepy::py {

namespace py = pybind11;

// CORE
void init_pyspline(py::module_&);

// Coordinate pointers
void init_coordinate_pointers(py::module_&);

// Knot Vector
void init_knot_vector(py::module_&);

// Parameter Space
void init_parameter_space(py::module_&);

// Extensions
void init_spline_extensions(py::module_& m);

// Reader
void init_spline_reader(py::module_&);

// Knot Insertion Matrix
void init_knot_insertion_matrix(py::module_&);

// Exporter
void init_spline_exporter(py::module_&);

// multipatch
void init_multipatch(py::module_& m);

} // namespace splinepy::py

namespace py = pybind11;

PYBIND11_MODULE(splinepy_core, m) {
  splinepy::py::init_pyspline(m);
  splinepy::py::init_coordinate_pointers(m);
  splinepy::py::init_knot_vector(m);
  splinepy::py::init_parameter_space(m);
  splinepy::py::init_spline_extensions(m);
  splinepy::py::init_spline_reader(m);
  splinepy::py::init_spline_exporter(m);
  splinepy::py::init_knot_insertion_matrix(m);
  splinepy::py::init_multipatch(m);

  // add some build configuration info
  m.def("build_type", []() {
#ifndef NDEBUG
    return "debug";
#else
  return "release";
#endif
  });
  m.def("is_minimal", []() {
#ifdef SPLINEPY_MORE
    return false;
#else
  return true;
#endif
  });
}
