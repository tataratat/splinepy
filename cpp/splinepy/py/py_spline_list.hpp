#pragma once

#include <cstdlib>

// pybind
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

#include <splinepy/py/py_spline.hpp>
#include <splinepy/utils/nthreads.hpp>

// PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<splinepy::py::PySpline>>);

namespace splinepy::py {

namespace py = pybind11;

using PySplineList = std::vector<std::shared_ptr<PySpline>>;

// evaluates splines of same para_dim and dim.
inline py::array_t<double> Evaluate(const PySplineList& splist,
                                    const py::array_t<double>& queries,
                                    const int nthreads) {
  const int& para_dim = splist[0]->para_dim_;
  const int& dim = splist[0]->dim_;

  CheckPyArrayShape(queries, {-1, para_dim}, true);

  // prepare input and output
  double* queries_ptr = static_cast<double*>(queries.request().ptr);
  const int n_splines = splist.size();
  const int n_queries = queries.size();
  const int n_total = n_splines * n_queries;
  py::array_t<double> evaluated({n_total, dim});
  double* evaluated_ptr = static_cast<double*>(evaluated.request().ptr);

  // each thread evaluates similar amount of queries from each spline
  auto evaluate = [&](int begin, int end) {
    const int n_common = std::div(end - begin, n_splines).quot; // floor
    const int spline_start = begin % n_splines;
    const int spline_end = (end - 1) % n_splines;
    const int query_offset = std::div(begin, n_splines).quot;

    // loop splines
    for (int i{}; i < n_splines; ++i) {
      const auto& core = *splist[i]->Core();
      const int output_offset = i * n_queries * dim;
      int query_start = query_offset;

      int n_additional{};
      if (i < spline_start) {
        n_additional -= 1;
        query_start += 1;
      }
      if (i <= spline_end) {
        n_additional += 1;
      }

      for (int j{query_start}; j < (query_start + n_common + n_additional);
           ++j) {
        core.Evaluate(&queries_ptr[j * para_dim],
                      &evaluated_ptr[output_offset + (j * dim)]);
      }
    }
  }
}

inline py::array_t<double> Derivative(const PySplineList& splist,
                                      const py::array_t<double>& queries,
                                      const py::array_t<int>& orders,
                                      const int nthreads) {}

inline py::array_t<double> Sample(const PySplineList& splist,
                                  const py::array_t<int> resolutions,
                                  const int nthreads) {}

inline std::shared_ptr<PySplineList>
ExtractBoundarySplines(const PySplineList& splist, const int nthreads) {}

inline py::array_t<double> BoundaryCenters(const PySplineList& splist,
                                           const int nthreads) {}

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
