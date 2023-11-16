#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <BSplineLib/ParameterSpaces/knot_vector.hpp>

#include "splinepy/utils/print.hpp"

namespace splinepy::py {

namespace py = pybind11;

void init_knot_vector(py::module_& m) {
  using KnotVector = bsplinelib::parameter_spaces::KnotVector;
  using IterType = typename KnotVector::Knots_::iterator;
  using KnotType = typename KnotVector::Knots_::value_type;
  using DiffType = typename KnotVector::Knots_::difference_type;
  using SizeType = typename KnotVector::Knots_::size_type;

  // implementations adapted from pybind11/stl_bind.h

  /// as python supports negative ids, this func brings negative ids to pos
  auto wrap_id = [](DiffType i, const SizeType n) {
    if (i < 0) {
      i += n;
    }
    if (i < 0 || (SizeType) i >= n) {
      throw py::index_error();
    }
    return i;
  };

  py::class_<KnotVector, std::shared_ptr<KnotVector>> klasse(m, "KnotVector");
  klasse
      .def(
          "__len__",
          [](const KnotVector& kv) { return kv.GetSize(); },
          "Returns size if len(knot_vector) is called.")
      .def(
          "__iter__",
          [](KnotVector& kv) {
            return py::make_iterator<
                py::return_value_policy::reference_internal,
                IterType,
                IterType,
                KnotType&>(kv.GetKnots().begin(), kv.GetKnots().end());
          },
          py::keep_alive<0, 1>(),
          "Support iterations of element references.")
      .def(
          "__getitem__",
          [wrap_id](KnotVector& kv, DiffType i) -> KnotType& {
            i = wrap_id(i, kv.GetKnots().size());
            return kv.GetKnots()[static_cast<SizeType>(i)];
          },
          py::return_value_policy::reference_internal,
          "int based __getitem__, which returns reference..")
      .def(
          "__getitem__",
          [](const KnotVector& kv, const py::slice& slice) -> py::list {
            std::size_t start{}, stop{}, step{}, slicelength{};
            const auto kv_size = kv.GetKnots().size();
            if (!slice.compute(kv_size, &start, &stop, &step, &slicelength)) {
              throw py::error_already_set();
            }
            py::list items{};
            for (std::size_t i{}; i < slicelength; ++i) {
              items.append(kv[start]);
              start += step;
            }
            return items;
          },
          py::arg("slice"),
          "support slice based __getitem__.")
      .def(
          "__setitem__",
          [wrap_id](KnotVector& kv, DiffType i, const KnotType knot) {
            i = wrap_id(i, kv.GetKnots().size());
            kv.UpdateKnot(i, knot);
          },
          "Single knot assignment / modification.")
      .def(
          "__setitem__",
          [](KnotVector& kv, const py::slice& slice, const py::list& value) {
            std::size_t start{}, stop{}, step{}, slicelength{};
            const auto kv_size = kv.GetKnots().size();
            if (!slice.compute(kv_size, &start, &stop, &step, &slicelength)) {
              throw py::error_already_set();
            }
            if (slicelength != value.size()) {
              splinepy::utils::PrintAndThrowError(
                  "Left and right hand size of slice assignment have "
                  "different sizes.");
            }
            auto& knots = kv.GetKnots();
            for (std::size_t i{}; i < slicelength; ++i) {
              knots[start] = py::cast<KnotType>(value[i]);
              start += step;
            }
            kv.ThrowIfTooSmallOrNotNonDecreasing();
          },
          "Multiple slice based element assignment.")
      .def("__repr__",
           [](const KnotVector& kv) {
             std::ostringstream s;
             s << "KnotVector [";
             int i{}, j{kv.GetSize() - 1};
             for (const auto& k : kv.GetKnots()) {
               s << k;
               if (i != j) {
                 s << ", ";
               }
               ++i;
             }
             s << "]";
             return s.str();
           })
      .def("scale",
           &KnotVector::Scale,
           py::arg("min"),
           py::arg("max"),
           "Scales knot vector with given [min, max].")
      .def("find_span",
           &KnotVector::FindSpan_,
           "Finds knot span of given parametric coordinate.")
      .def(
          "numpy",
          [](const KnotVector& kv) {
            py::array_t<KnotType> arr(kv.GetSize());
            KnotType* arr_ptr = static_cast<KnotType*>(arr.request().ptr);
            for (int i{}; i < kv.GetSize(); ++i) {
              arr_ptr[i] = kv[i];
            }
            return arr;
          },
          "Returns copy of knot vectors as numpy array.")
      .def(
          "__array__",
          [](const KnotVector& kv, [[maybe_unused]] py::args dtype_ignored) {
            py::array_t<KnotType> arr(kv.GetSize());
            KnotType* arr_ptr = static_cast<KnotType*>(arr.request().ptr);
            for (int i{}; i < kv.GetSize(); ++i) {
              arr_ptr[i] = kv[i];
            }
            return arr;
          },
          "Similar to numpy(), but supports np.array or np.asarray based "
          "creation.")
      .def(
          "unique",
          [](const KnotVector& kv) -> py::array_t<KnotType> {
            const KnotVector::Knots_& uniques = kv.GetUniqueKnots();
            py::array_t<KnotType> arr(uniques.size());
            KnotType* arr_ptr = static_cast<KnotType*>(arr.request().ptr);
            for (int i{}; i < static_cast<int>(uniques.size()); ++i) {
              arr_ptr[i] = uniques[i];
            }
            return arr;
          },
          "Returns multiplicities of unique knots")
      .def(
          "multiplicities",
          [](const KnotVector& kv) -> py::array_t<int> {
            const bsplinelib::Vector<int> multiplicities =
                kv.DetermineMultiplicities();
            py::array_t<int> arr(multiplicities.size());
            int* arr_ptr = static_cast<int*>(arr.request().ptr);
            for (int i{}; i < static_cast<int>(multiplicities.size()); ++i) {
              arr_ptr[i] = multiplicities[i];
            }
            return arr;
          },
          "Returns multiplicities of unique knots");
}

} // namespace splinepy::py
