#include <pybind11/pybind11.h>

// core_spline
namespace splinepy::py::init {

namespace py = pybind11;

// CORE
void init_core_spline(py::module_&);

// Coordinate pointers
void init_coordinate_pointers(py::module_&);

// Knot Vector
void init_knot_vector(py::module_&);

// Extensions
void init_spline_extensions(py::module_& m);

// Reader
void init_reader(py::module_&);

// Knot Insertion Matrix
void init_knot_insertion_matrix(py::module_&);

// Exporter
void init_exporter(py::module_&);

// fitting
void init_fitting(py::module_&);

// multi_patch
void init_multi_patch(py::module_& m);

} // namespace splinepy::py::init

namespace py = pybind11;

PYBIND11_MODULE(splinepy_core, m) {
  splinepy::py::init::init_core_spline(m);
  splinepy::py::init::init_coordinate_pointers(m);
  splinepy::py::init::init_knot_vector(m);
  splinepy::py::init::init_spline_extensions(m);
  splinepy::py::init::init_reader(m);
  splinepy::py::init::init_exporter(m);
  splinepy::py::init::init_fitting(m);
  splinepy::py::init::init_knot_insertion_matrix(m);
  splinepy::py::init::init_multi_patch(m);
}
