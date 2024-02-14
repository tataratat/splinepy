#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <BSplineLib/ParameterSpaces/knot_vector.hpp>
#include <BSplineLib/ParameterSpaces/parameter_space.hpp>

#include "splinepy/utils/print.hpp"

namespace splinepy::py {

namespace py = pybind11;

void init_parameter_space(py::module_& m) {
  using PSpace = bsplinelib::parameter_spaces::ParameterSpaceBase;
  using KV = bsplinelib::parameter_spaces::KnotVector;
  using KVPtr = std::shared_ptr<KV>;

  // implementations adapted from pybind11/stl_bind.h
  /// as python supports negative ids, this func brings negative ids to pos
  auto wrap_id = [](int i, int n) -> int {
    if (i < 0) {
      i += n;
    }
    if (i < 0 || i >= n) {
      throw py::index_error();
    }
    return i;
  };

  auto to_list = [](const PSpace& p) -> py::list {
    // if this were to be a list, this would mean a shallow copy.
    // however, within splinepy's implementation, copy() means deepcopy.
    // so, deepcopy it is.
    py::list kvs(p.ParaDim());
    for (int i{}; i < p.ParaDim(); ++i) {
      KVPtr kv = p.GetKnotVector(i);
      py::array_t<double> out_kv(kv->GetSize());
      std::copy_n(kv->GetKnots().data(),
                  kv->GetKnots().size(),
                  static_cast<double*>(out_kv.request().ptr));
      kvs[i] = out_kv;
    }
    return kvs;
  };

  py::class_<PSpace, std::shared_ptr<PSpace>> klasse(m, "ParameterSpace");

  klasse
      .def(
          "__len__",
          [](const PSpace& p) { return p.ParaDim(); },
          "Returns number of parameter dimension.")

      .def(
          "__iter__",
          [](PSpace& p) {
            return py::make_iterator<
                py::return_value_policy::reference_internal,
                KVPtr*,
                KVPtr*,
                KVPtr&>(p.KnotVectorsBegin(), p.KnotVectorsEnd());
          },
          py::keep_alive<0, 1>(),
          "Support iterations of element references.")

      .def(
          "__getitem__",
          [wrap_id](PSpace& p, int i) -> KVPtr& {
            i = wrap_id(i, p.ParaDim());
            return p.GetKnotVector(i);
          },
          py::return_value_policy::reference_internal,
          "int based __getitem__, which returns reference.")
      .def(
          "__getitem__",
          [](const PSpace& p, const py::slice& slice) -> py::list {
            std::size_t start{}, stop{}, step{}, slicelength{};
            if (!slice.compute(static_cast<std::size_t>(p.ParaDim()),
                               &start,
                               &stop,
                               &step,
                               &slicelength)) {
              throw py::error_already_set();
            }
            py::list items{};
            for (std::size_t i{}; i < slicelength; ++i) {
              items.append(p.GetKnotVector(static_cast<int>(start)));
              start += step;
            }
            return items;
          },
          py::arg("slice"),
          "support slice based __getitem__.")
      .def(
          "__setitem__",
          [wrap_id](PSpace& p, int i, const KVPtr new_kv) {
            i = wrap_id(i, p.ParaDim());
            auto& p_kv = p.GetKnotVector(i);
            if (p_kv->GetSize() != new_kv->GetSize()) {
              splinepy::utils::PrintAndThrowError(
                  "Size mismatch of lhs & rhs knot vectors");
            }
            // update is simple assignment - note that this is reference
            // this is a same behavir as list
            p_kv = new_kv;
          },
          "Single knot vector assignment with another knot vector.")
      .def(
          "__setitem__",
          [wrap_id](PSpace& p, int i, const py::array_t<double>& new_kv) {
            i = wrap_id(i, p.ParaDim());
            auto& p_kv = p.GetKnotVector(i);
            const int kv_size = p_kv->GetSize();
            if (kv_size != new_kv.size()) {
              splinepy::utils::PrintAndThrowError(
                  "Size mismatch of lhs & rhs knot vectors");
            }
            // update is value copy
            double* p_data = p_kv->GetKnots().data();
            std::copy_n(static_cast<double*>(new_kv.request().ptr),
                        kv_size,
                        p_data);

            // sanity check
            p_kv->ThrowIfTooSmallOrNotNonDecreasing();
          },
          "Single knot vector assignment with an array")
      .def("__add__",
           [to_list](const PSpace& p, py::list& next) {
             return to_list(p) + next;
           })
      .def("__radd__",
           [to_list](const PSpace& p, py::list& next) {
             return next + to_list(p);
           })
      .def("copy", [to_list](const PSpace& p) { return to_list(p); });
}

} // namespace splinepy::py
