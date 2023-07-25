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

// fitting
void init_fitting(py::module_&);

// multi_patch
// void init_multi_patch(py::module_& m);
void add_multi_patch(py::module_& m);

} // namespace splinepy::py

namespace py = pybind11;

PYBIND11_MODULE(splinepy_core, m) {
  splinepy::py::init_pyspline(m);
  splinepy::py::init_coordinate_pointers(m);
  splinepy::py::init_knot_vector(m);
  splinepy::py::init_spline_extensions(m);
  splinepy::py::init_spline_reader(m);
  splinepy::py::init_spline_exporter(m);
  splinepy::py::init_fitting(m);
  splinepy::py::init_knot_insertion_matrix(m);
  // splinepy::py::init_multi_patch(m);
  splinepy::py::add_multi_patch(m);
}
