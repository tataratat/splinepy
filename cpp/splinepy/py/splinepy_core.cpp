#include <pybind11/pybind11.h>

namespace py = pybind11;

// core_spline
void init_core_spline(py::module_&);

// Reader
void init_reader(py::module_&);

// Exporter
void init_exporter(py::module_&);

PYBIND11_MODULE(splinepy_core, m) {

  init_core_spline(m);
  init_reader(m);
  init_exporter(m);
}
