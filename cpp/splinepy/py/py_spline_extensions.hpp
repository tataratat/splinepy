/// collection of functions that take spline as parameter
/// but either not applicable for all splines or only applicable to
/// specific use cases.
#pragma once

// pybind11
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <splinepy/py/py_spline.hpp>

namespace splinepy::py {

namespace py = pybind11;

/// (multiple) knot insertion, single dimension
inline void InsertKnots(std::shared_ptr<PySpline>& spline,
                        int para_dim,
                        py::array_t<double> knots) {
  double* knots_ptr = static_cast<double*>(knots.request().ptr);
  const int n_request = knots.size();

  for (int i{}; i < n_request; ++i) {
    spline->Core()->SplinepyInsertKnot(para_dim, knots_ptr[i]);
  }
}

/// (multiple) knot removal, single dimension
inline py::list RemoveKnots(std::shared_ptr<PySpline>& spline,
                            int para_dim,
                            py::array_t<double> knots,
                            double tolerance) {
  double* knots_ptr = static_cast<double*>(knots.request().ptr);
  const int n_request = knots.size();

  py::list successful;
  for (int i{}; i < n_request; ++i) {
    successful.append(
        spline->Core()->SplinepyRemoveKnot(para_dim, knots_ptr[i], tolerance));
  }

  return successful;
}

/// spline multiplication - currently only for bezier
inline std::shared_ptr<PySpline> Multiply(const std::shared_ptr<PySpline>& a,
                                          const std::shared_ptr<PySpline>& b) {
  // performs runtime checks and throws error
  return std::make_shared<PySpline>(a->Core()->SplinepyMultiply(b->Core()));
}

/// spline addition - currently only for bezier
inline std::shared_ptr<PySpline> Add(const std::shared_ptr<PySpline>& a,
                                     const std::shared_ptr<PySpline>& b) {
  // performs runtime checks and throws error
  return std::make_shared<PySpline>(a->Core()->SplinepyAdd(b->Core()));
}

/// spline composition - currently only for bezier
inline std::shared_ptr<PySpline>
Compose(const std::shared_ptr<PySpline>& inner,
        const std::shared_ptr<PySpline>& outer) {
  // performs runtime checks and throws error
  return std::make_shared<PySpline>(
      inner->Core()->SplinepyCompose(outer->Core()));
}

/// spline derivative spline
inline std::shared_ptr<PySpline>
DerivativeSpline(const std::shared_ptr<PySpline>& spline,
                 py::array_t<int> orders) {
  CheckPyArraySize(orders, spline->para_dim_);

  int* orders_ptr = static_cast<int*>(orders.request().ptr);
  return std::make_shared<PySpline>(
      spline->Core()->SplinepyDerivativeSpline(orders_ptr));
}

/// spline split - returns py::list of PySplines
inline py::list Split(const std::shared_ptr<PySpline>& spline,
                      int p_dim,
                      py::array_t<double> locations) {
  // make sure they are sorted
  std::vector<double> locs(locations.size());
  std::copy_n(static_cast<double*>(locations.request().ptr),
              locations.size(),
              locs.begin());
  // sort
  std::sort(locs.begin(), locs.end());

  // split and append
  py::list splitted;
  // very first
  auto tmp_splitted = spline->Core()->SplinepySplit(p_dim, locs[0]);
  splitted.append(std::make_shared<PySpline>(tmp_splitted[0]));
  for (std::size_t i{1}; i < locs.size(); ++i) {
    const double locs_with_offset =
        (locs[i] - locs[i - 1]) / (1. - locs[i - 1]);
    tmp_splitted = tmp_splitted[1]->SplinepySplit(p_dim, locs_with_offset);
    splitted.append(std::make_shared<PySpline>(tmp_splitted[0]));
  }
  // very last
  splitted.append(std::make_shared<PySpline>(tmp_splitted[1]));

  return splitted;
}

/// bezier patch extraction.
inline py::list ExtractBezierPatches(const std::shared_ptr<PySpline>& spline) {
  const auto sp_patches = spline->Core()->SplinepyExtractBezierPatches();
  py::list patches;
  for (const auto& p : sp_patches) {
    patches.append(std::make_shared<PySpline>(p));
  }
  return patches;
}

/// boundary spline extraction
inline py::list ExtractBoundaries(const std::shared_ptr<PySpline>& spline,
                                  const py::array_t<int>& boundary_ids) {
  // Init return value
  py::list boundary_splines{};
  const int n_boundaries = boundary_ids.size();
  int* bid_ptr = static_cast<int*>(boundary_ids.request().ptr);
  if (boundary_ids.size() == 0) {
    for (int i{}; i < spline->para_dim_ * 2; ++i) {
      boundary_splines.append(std::make_shared<PySpline>(
          spline->Core()->SplinepyExtractBoundary(i)));
    }
  } else {
    const int max_bid = spline->para_dim_ * 2 - 1;
    for (int i{}; i < n_boundaries; ++i) {
      const int& bid = bid_ptr[i];
      if (bid < 0 || bid > max_bid) {
        splinepy::utils::PrintAndThrowError("Requested Boundary ID :",
                                            bid,
                                            "exceeds admissible range.");
      }
      boundary_splines.append(std::make_shared<PySpline>(
          spline->Core()->SplinepyExtractBoundary(bid_ptr[i])));
    }
  }

  return boundary_splines;
}

/// extract a single physical dimension from a spline
inline std::shared_ptr<PySpline>
ExtractDim(const std::shared_ptr<PySpline>& spline, int phys_dim) {
  return std::make_shared<PySpline>(
      spline->Core()->SplinepyExtractDim(phys_dim));
}

/// composition derivative
inline std::shared_ptr<PySpline>
CompositionDerivative(const std::shared_ptr<PySpline>& outer,
                      const std::shared_ptr<PySpline>& inner,
                      const std::shared_ptr<PySpline>& inner_derivative) {
  return std::make_shared<PySpline>(
      outer->Core()->SplinepyCompositionDerivative(inner->Core(),
                                                   inner_derivative->Core()));
}

/// returns a spline with knot vectors.
/// if the spline already has knots, it returns the same spline
/// else, returns a same spline with knots
inline std::shared_ptr<PySpline>
SameSplineWithKnotVectors(std::shared_ptr<PySpline>& spline) {
  // early exit if the spline has knot vectors already
  if (spline->HasKnotVectors()) {
    return spline;
  }

  py::dict props = spline->CurrentCoreProperties();

  // based on degrees, generate knot vectors
  py::array_t<int> degrees = py::cast<py::array_t<int>>(props["degrees"]);
  int* degrees_ptr = static_cast<int*>(degrees.request().ptr);
  py::list kvs;
  for (int i{}; i < degrees.size(); ++i) {
    // prepare kv and get required number of repeating knots
    const int n_knot_repeat = degrees_ptr[i] + 1;
    py::array_t<double> kv(n_knot_repeat * 2);
    double* kv_ptr = static_cast<double*>(kv.request().ptr);

    // defined knots with 0 and 1
    for (int j{}; j < n_knot_repeat; ++j) {
      kv_ptr[j] = 0.;
      kv_ptr[n_knot_repeat + j] = 1.;
    }

    kvs.append(kv);
  }

  // update knot vectors to dict spline
  props["knot_vectors"] = kvs;

  return std::make_shared<PySpline>(props);
}

/// @brief Evaluate Splines at boundary face centers
/// @return numpy array with results
inline py::array_t<double>
EvaluateBoundaryCenters(std::shared_ptr<PySpline>& spline) {
  // prepare output
  const int& para_dim_ = spline->para_dim_;
  const int& dim_ = spline->dim_;
  const int n_faces = 2 * para_dim_;
  std::vector<double> queries(n_faces * para_dim_);
  py::array_t<double> face_centers(n_faces * dim_);
  double* face_centers_ptr = static_cast<double*>(face_centers.request().ptr);

  // Get parametric bounds
  // They are given back in the order [min_0, min_1,...,max_0, max_1...]
  std::vector<double> bounds(para_dim_ * 2);
  spline->Core()->SplinepyParametricBounds(bounds.data());

  // Set parametric coordinates
  for (int i{}; i < para_dim_; i++) {
    for (int j{}; j < para_dim_; j++) {
      if (i == j) {
        queries[2 * i * para_dim_ + j] = bounds[j];
        queries[(2 * i + 1) * para_dim_ + j] = bounds[j + para_dim_];
      } else {
        queries[2 * i * para_dim_ + j] =
            .5 * (bounds[j] + bounds[j + para_dim_]);
        queries[(2 * i + 1) * para_dim_ + j] = queries[2 * i * para_dim_ + j];
      }
    }
  }
  for (int i{}; i < n_faces; ++i) {
    spline->Core()->SplinepyEvaluate(&queries.data()[i * para_dim_],
                                     &face_centers_ptr[i * dim_]);
  }

  face_centers.resize({n_faces, dim_});
  return face_centers;
}

/// returns core spline's ptr address
inline intptr_t CoreId(const std::shared_ptr<PySpline>& spline) {
  return reinterpret_cast<intptr_t>(spline->Core().get());
}

/// reference count of core spline
inline int CoreRefCount(const std::shared_ptr<PySpline>& spline) {
  return spline->Core().use_count();
}

/// have core? A non error raising checker
inline bool HaveCore(const std::shared_ptr<PySpline>& spline) {
  return (spline->c_spline_) ? true : false;
}

/// Overwrite core with a nullptr and assign neg values to dims
inline void AnnulCore(std::shared_ptr<PySpline>& spline) {
  spline->c_spline_ = nullptr;
  spline->para_dim_ = -1;
  spline->dim_ = -1;
}

inline void add_spline_extensions(py::module& m) {
  m.def("insert_knots",
        &splinepy::py::InsertKnots,
        py::arg("spline"),
        py::arg("para_dim"),
        py::arg("knots"));
  m.def("remove_knots",
        &splinepy::py::RemoveKnots,
        py::arg("spline"),
        py::arg("para_dim"),
        py::arg("knots"),
        py::arg("tolerance"));
  m.def("multiply", &splinepy::py::Multiply, py::arg("a"), py::arg("b"));
  m.def("add", &splinepy::py::Add, py::arg("a"), py::arg("b"));
  m.def("compose", &splinepy::py::Compose, py::arg("outer"), py::arg("inner"));
  m.def("derivative_spline",
        &splinepy::py::DerivativeSpline,
        py::arg("spline"),
        py::arg("orders"));
  m.def("split",
        &splinepy::py::Split,
        py::arg("spline"),
        py::arg("para_dim"),
        py::arg("locations"));
  m.def("extract_bezier_patches",
        &splinepy::py::ExtractBezierPatches,
        py::arg("spline"));
  m.def("extract_boundaries",
        &splinepy::py::ExtractBoundaries,
        py::arg("spline"),
        py::arg("boundary_ids"));
  m.def("extract_dim",
        &splinepy::py::ExtractDim,
        py::arg("spline"),
        py::arg("phys_dim"));
  m.def("composition_derivative",
        &splinepy::py::CompositionDerivative,
        py::arg("outer"),
        py::arg("inner"),
        py::arg("inner_derivative"));
  m.def("same_spline_with_knot_vectors",
        &splinepy::py::SameSplineWithKnotVectors,
        py::arg("spline"));
  m.def("boundary_centers",
        &splinepy::py::EvaluateBoundaryCenters,
        py::arg("spline"));
  m.def("core_id", &splinepy::py::CoreId, py::arg("spline"));
  m.def("core_ref_count", &splinepy::py::CoreRefCount, py::arg("spline"));
  m.def("have_core", &splinepy::py::HaveCore, py::arg("spline"));
  m.def("annul_core", &splinepy::py::AnnulCore, py::arg("spline"));
  ;
}

} // namespace splinepy::py
