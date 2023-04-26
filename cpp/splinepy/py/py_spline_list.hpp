#pragma once

// pybind
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

#include <splinepy/py/py_spline.hpp>

PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<splinepy::py::PySpline>>);

namespace splinepy::py {

namespace py = pybind11;

using PySplineList = std::vector<std::shared_ptr<PySpline>>;


inline py::array_t<double> Evaluate(const PySplineList& splist, const int nthreads) {
}

inline py::array_t<double> EvaluateUsingSameBasis(const PySplineList& splist, const int nthreads) {
}

inline py::array_t<double> Derivative(const PySplineList& splist, const int nthreads) {
}

inline py::array_t<double> DerivativeUsingSameBasis(const PySplineList& splist, const int nthreads) {
}

inline py::array_t<double> Sample(const PySplineList& splist, const int nthreads) {
}



/// bind vector of PySpline and add some deprecated cpp functions that maybe
/// nice to have
inline void add_spline_list_pyclass(py::module& m) {

  // use shared_ptr as holder
  py::bind_vector<PySplineList, std::shared_ptr<PySplineList>>(m, "SplineList");

  m.def("reserve_list",
        [](PySplineList& splist, py::ssize_t size) { splist.reserve(size); });
  m.def("resize_list",
        [](PySplineList& splist, py::ssize_t size) { splist.resize(size); });
  m.def("shrink_to_fit_list",
        [](PySplineList& splist) { splist.shrink_to_fit(); });
  m.def("swap_list", [](PySplineList& splist_a, PySplineList& splist_b) {
    splist_a.swap(splist_b);
  });
}

} // namespace splinepy::py
