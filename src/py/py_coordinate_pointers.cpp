#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "splinepy/utils/coordinate_pointers.hpp"
#include "splinepy/utils/print.hpp"

namespace splinepy::py {

namespace py = pybind11;

/// Wraps ControlPointPointers and WeightPointers for python.
/// Before passing syncing values, this wrappers perform size checks.

template<typename PointersType>
int Len(const PointersType& pointers) {
  return pointers.Len();
}

template<typename PointersType>
int Dim(const PointersType& pointers) {
  return pointers.Dim();
}

template<typename PointersType, bool same_sized_values>
void SetRows(PointersType& pointers,
             const py::array_t<int> ids,
             const py::array_t<double>& values) {
  const int n_rows = ids.size();

  if constexpr (same_sized_values) {
    if (pointers.Len() * pointers.Dim() != values.size()) {
      splinepy::utils::PrintAndThrowError("Size mismatch. Expecting",
                                          pointers.Len() * pointers.Dim(),
                                          "but, values are ",
                                          values.size());
    }
  } else {
    if (n_rows * pointers.Dim() != values.size()) {
      splinepy::utils::PrintAndThrowError("Size mismatch. Expecting",
                                          pointers.Dim() * ids.size(),
                                          "but, values are ",
                                          values.size());
    }
  }

  if (pointers.Dim() != values.shape(1)) {
    splinepy::utils::PrintAndThrowError("Dimension mismatch. Expecting",
                                        pointers.Dim(),
                                        "but values are",
                                        values.shape(1));
  }
  pointers.template SetRows<same_sized_values>(
      static_cast<int*>(ids.request().ptr),
      n_rows,
      static_cast<double*>(values.request().ptr));
}

template<typename PointersType>
void Sync(PointersType& pointers, const py::array_t<double>& values) {
  if (pointers.Len() * pointers.Dim() != values.size()) {
    splinepy::utils::PrintAndThrowError("Size mismatch, Expecting",
                                        pointers.Len() * pointers.Dim(),
                                        "but values are",
                                        values.size());
  }
  pointers.Sync(static_cast<double*>(values.request().ptr));
}

template<typename PointersType>
std::shared_ptr<PointersType> SubSet(PointersType& pointers,
                                     const py::array_t<int>& ids) {
  return pointers.SubSet(static_cast<int*>(ids.request().ptr), ids.size());
}

void init_coordinate_pointers(py::module& m) {
  using ControlPointPointers = splinepy::utils::ControlPointPointers;
  using WeightPointers = splinepy::utils::WeightPointers;

  py::class_<ControlPointPointers, std::shared_ptr<ControlPointPointers>>
      klasse(m, "ControlPointPointers");
  klasse.def("len", &Len<ControlPointPointers>)
      .def("dim", &Dim<ControlPointPointers>)
      .def(
          "set_row",
          [](ControlPointPointers& pointers,
             const int id,
             const py::array_t<double>& values) {
            if (pointers.Dim() != values.size()) {
              splinepy::utils::PrintAndThrowError("Size mismatch. Expecting",
                                                  pointers.dim_,
                                                  "but, values are ",
                                                  values.size(),
                                                  ".");
            }
            pointers.SetRow(id, static_cast<double*>(values.request().ptr));
          },
          py::arg("id"),
          py::arg("values"))
      .def("set_rows",
           &SetRows<ControlPointPointers, false>,
           py::arg("id"),
           py::arg("values"))
      .def("sync_rows",
           &SetRows<ControlPointPointers, true>,
           py::arg("id"),
           py::arg("values"))
      .def("sync", &Sync<ControlPointPointers>, py::arg("values"));

  py::class_<WeightPointers, std::shared_ptr<WeightPointers>> klasse_weight(
      m,
      "WeightPointers");
  klasse_weight.def("len", &Len<WeightPointers>)
      .def("dim", &Dim<WeightPointers>)
      .def("set_row", &WeightPointers::SetRow, py::arg("id"), py::arg("values"))
      .def("set_rows",
           &SetRows<WeightPointers, false>,
           py::arg("id"),
           py::arg("values"))
      .def("sync_rows",
           &SetRows<WeightPointers, true>,
           py::arg("id"),
           py::arg("values"))
      .def("sync", &Sync<WeightPointers>, py::arg("values"));
}

} // namespace splinepy::py
