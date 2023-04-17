#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <splinepy/splines/splinepy_base.hpp>
#include <splinepy/utils/print.hpp>
#include <splinepy/utils/reference.hpp>

namespace splinepy::py {

namespace py = pybind11;

/// Access reference of each control point coefficients as if they were
/// layed out in contiguous manner.
/// excess by global index.
/// probably useful to wrap this once with numpy __getitem__
/// we name this coordinates, since rational splines will have weighted
/// control points.
using CoordinateReferences =
    splinepy::splines::SplinepyBase::CoordinateReferences_;

inline void add_coordinate_references_pyclass(py::module_& m) {
  py::class_<CoordinateReferences, std::shared_ptr<CoordinateReferences>>
      klasse(m, "CoordinateReferences");
  klasse.def(py::init<>())
      .def("__len__", [](const CoordinateReferences& r) { return r.size(); })
      .def("__setitem__",
           [](CoordinateReferences& r,
              const py::array_t<int> ids,
              const py::array_t<double> values) {
             // size match check
             const py::ssize_t i_size = ids.size();
             const py::ssize_t v_size = values.size();

             if (i_size != v_size) {
               splinepy::utils::PrintAndThrowError(
                   "array size mismatch between ids (",
                   i_size,
                   ")  and values (",
                   v_size,
                   ")");
             }
             // get ptr
             int* ids_ptr = static_cast<int*>(ids.request().ptr);
             double* values_ptr = static_cast<double*>(values.request().ptr);

             // assign
             for (py::ssize_t i{}; i < i_size; ++i) {
               r[ids_ptr[i]].value_ = values_ptr[i];
             }
           })
      .def("broadcast_scalar",
           [](CoordinateReferences& r,
              const py::array_t<int> ids,
              const double value) {
             const py::ssize_t i_size = ids.size();
             int* ids_ptr = static_cast<int*>(ids.request().ptr);
             for (py::ssize_t i{}; i < i_size; ++i) {
               r[ids_ptr[i]].value_ = value;
             }
           })
      .def("numpy", [](const CoordinateReferences& r) {
        const int r_size = r.size();
        py::array_t<double> copied(r_size);
        double* copied_ptr = static_cast<double*>(copied.request().ptr);
        for (int i{}; i < r_size; ++i) {
          copied_ptr[i] = r[i].value_;
        }
        return copied;
      });
}

} // namespace splinepy::py
