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
