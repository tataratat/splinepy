#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "third_party/uff/src/uff.hpp"

namespace splinepy::py {

namespace py = pybind11;

/**
 * @brief
 *
 * @param points
 * @param tolerance
 * @param stable
 * @return py::tuple
 */
py::tuple uffpy(py::array_t<double> points, double tolerance, bool stable) {
  // Access points
  double* p_buf_ptr = static_cast<double*>(points.request().ptr);
  int npoints = points.request().shape[0];
  int pdim = points.request().shape[1];

  // selecting metric for you.
  std::vector<double> metric(pdim, 1.);

  // prepare output arrays
  py::array_t<int> np_newpointmasks(npoints);
  py::buffer_info npm_buf = np_newpointmasks.request();
  int* newpointmasks = static_cast<int*>(npm_buf.ptr);

  py::array_t<int> np_inverse(npoints);
  py::buffer_info i_buf = np_inverse.request();
  int* inverse = static_cast<int*>(i_buf.ptr);

  // prepare additional input vars
  int nnewpoints;

  // prepare temp array to store newpoints.
  // we call this thing phil
  py::array_t<double> np_newpoints({nnewpoints, pdim});

  // Dialogue
  // ---------
  // pyuff - (calls uff)
  // uff - Guten Tag uff am Apparat
  // pyuff - einmal uff bitte
  // uff - ::uff!
  uff::uff(p_buf_ptr,
           npoints,
           pdim,
           metric.data(),
           tolerance,
           stable,
           static_cast<double*>(np_newpoints.request().ptr),
           newpointmasks,
           nnewpoints,
           inverse);

  // real newpoints
  np_newpoints.resize({nnewpoints, pdim}, false);

  return py::make_tuple(np_newpoints, np_newpointmasks, np_inverse);
}

inline void add_uffpy(py::module& m) {
  m.def("unique_vertices_uff",
        &splinepy::py::uffpy,
        py::arg("queries"),
        py::arg("tolerance"),
        py::arg("stable") = true);
}

} // namespace splinepy::py
