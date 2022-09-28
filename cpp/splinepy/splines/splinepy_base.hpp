#pragma once

#include <memory>
#include <vector>

// pybind11
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <splinepy/utils/print.hpp>

namespace splinepy::splines {

namespace py = pybind11;

using namespace splinelib::sources;

/// Spline base to enable dynamic use of template splines.
/// Member functions are prepended with "RawPtr".
class SplinepyBase {
public:
  SplinepyBase() = default;

  static std::shared_ptr<SplinepyBase> Create(
      const int para_dim = 0,
      const int dim = 0,
      const double* degrees = nullptr,
      const std::vector<std::vector<double>>* knot_vectors = nullptr,
      const double* control_points = nullptr,
      const double* weights = nullptr) {
    if (!degrees || !control_points) {
      splinepy::utils::PrintAndThrowError(
          "Not Enough information to create any spline."
      );
    }

    if (!knot_vectors) {
      // at least we need degrees and cps.
      // @jzwar this would be a good place to check valid input

      if (!weights) {
        //return CreateBezier(degrees, control_points);
        return CreateBezier();
      } else {
        //return CreateRationalBezier(degrees, control_points, weights);
        return CreateRationalBezier();
      }
    } else {
      if (!weights) {
        //return CreateBSpline(degrees, knot_vectors, control_points);
        return CreateBSpline();
      } else {
        return CreateNurbs(para_dim,
                           dim,
                           degrees,
                           knot_vectors,
                           control_points,
                           weights);
      }
    }
  };

  // Implemented in splinepy/splines/bezier.hpp
  static std::shared_ptr<SplinepyBase> CreateBezier() {};
  // Implemented in splinepy/splines/rational_bezier.hpp
  static std::shared_ptr<SplinepyBase> CreateRationalBezier() {};
  // Implemented in splinepy/splines/b_spline.hpp
  static std::shared_ptr<SplinepyBase> CreateBSpline() {};
  /// Implemented in splinepy/splines/nurbs.hpp
  static std::shared_ptr<SplinepyBase> CreateNurbs(
      const int para_dim,
      const int dim,
      const double* degrees,
      const std::vector<std::vector<double>>* knot_vectors,
      const double* control_points,
      const double* weights);

  virtual constexpr int ParaDim() const = 0;
  virtual constexpr int Dim() const = 0;
  virtual std::string WhatAmI() const = 0;
  virtual int NumberOfControlPoints() const = 0;
  /// Extract core spline properties. Similar to previous update_p
  virtual void RawPtrCurrentProperties(
      double* degrees,
      std::vector<std::vector<double>>* knot_vectors,
      double* control_points,
      double* weights
  ) const = 0;

  virtual void RawPtrParametricBounds(double* p_bounds) const {
    splinepy::utils::PrintAndThrowError(
        "RawPtrParametricBounds not implemented for",
        WhatAmI()
    );

  };

  virtual void RawPtrEvaluate(double* para_coord, double* evaluated) {
    splinepy::utils::PrintAndThrowError(
        "RawPtrEvaluate not implemented for",
        WhatAmI()
    );
  };

  virtual void RawPtrDerivative(double* para_coord,
                                double* orders,
                                double* derived) {
    splinepy::utils::PrintAndThrowError(
        "RawPtrDerivative not implemented for",
        WhatAmI()
    );
  };

  virtual void RawPtrElevateDegrees(int* para_dims) {
    splinepy::utils::PrintAndThrowError(
        "RawPtrElevateDegrees not implemented for",
        WhatAmI()
    );
  };

  virtual void RawPtrReduceDegrees(int* para_dims, double tolerance) {
    splinepy::utils::PrintAndThrowError(
        "RawPtrReduceDegrees not implemented for",
        WhatAmI()
    );
  };

  virtual void RawPtrInsertKnots(int para_dim, int* knots) {
    splinepy::utils::PrintAndThrowError(
        "RawPtrInsertKnots not implemented for",
        WhatAmI()
    );
  };

  virtual void RawPtrRemoveKnots(int para_dim, int* knots) {
    splinepy::utils::PrintAndThrowError(
        "RawPtrRemoveKnots not implemented for",
        WhatAmI()
    );
  };
};

} // namespace splinepy::splines
