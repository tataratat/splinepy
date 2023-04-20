#include <pybind11/pybind11.h>

// core_spline
namespace splinepy::py::init {

namespace py = pybind11;

// Coordinate ref
void init_coordinate_references(py::module_&);

// CORE
void init_core_spline(py::module_&);

// Extensions
void init_spline_extensions(py::module_& m);

// Reader
void init_reader(py::module_&);

// Exporter
void init_exporter(py::module_&);

// fitting
void init_fitting(py::module_&);

} // namespace splinepy::py::init

namespace py = pybind11;

// PYBIND11_MAKE_OPAQUE(splinepy::py::CoordinateReferences);

PYBIND11_MODULE(splinepy_core, m) {
  splinepy::py::init::init_coordinate_references(m);
  splinepy::py::init::init_core_spline(m);
  splinepy::py::init::init_spline_extensions(m);
  splinepy::py::init::init_reader(m);
  splinepy::py::init::init_exporter(m);
  splinepy::py::init::init_fitting(m);
}
