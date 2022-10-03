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
/// Member functions are prepended with "Splinepy".
class SplinepyBase {
public:
  SplinepyBase() = default;

  static std::shared_ptr<SplinepyBase>
  SplinepyCreate(const int para_dim = 0,
                 const int dim = 0,
                 const double* degrees = nullptr,
                 const std::vector<std::vector<double>>* knot_vectors = nullptr,
                 const double* control_points = nullptr,
                 const double* weights = nullptr) {
    if (!degrees || !control_points) {
      splinepy::utils::PrintAndThrowError(
          "Not Enough information to create any spline.");
    }

    if (!knot_vectors) {
      // at least we need degrees and cps.
      // @jzwar this would be a good place to check valid input

      if (!weights) {
        return SplinepyCreateBezier(para_dim, dim, degrees, control_points);
      } else {
        return SplinepyCreateRationalBezier(para_dim,
                                            dim,
                                            degrees,
                                            control_points,
                                            weights);
      }
    } else {
      if (!weights) {
        return SplinepyCreateBSpline(para_dim,
                                     dim,
                                     degrees,
                                     knot_vectors,
                                     control_points);
      } else {
        return SplinepyCreateNurbs(para_dim,
                                   dim,
                                   degrees,
                                   knot_vectors,
                                   control_points,
                                   weights);
      }
    }
  };

  // Implemented in splinepy/splines/bezier.hpp
  static std::shared_ptr<SplinepyBase>
  SplinepyCreateBezier(const int para_dim,
                       const int dim,
                       const double* degrees,
                       const double* control_points);

  // Implemented in splinepy/splines/rational_bezier.hpp
  static std::shared_ptr<SplinepyBase>
  SplinepyCreateRationalBezier(const int para_dim,
                               const int dim,
                               const double* degrees,
                               const double* control_points,
                               const double* weights);
  // Implemented in splinepy/splines/b_spline.hpp
  static std::shared_ptr<SplinepyBase>
  SplinepyCreateBSpline(const int para_dim,
                        const int dim,
                        const double* degrees,
                        const std::vector<std::vector<double>>* knot_vectors,
                        const double* control_points);

  /// Implemented in splinepy/splines/nurbs.hpp
  static std::shared_ptr<SplinepyBase>
  SplinepyCreateNurbs(const int para_dim,
                      const int dim,
                      const double* degrees,
                      const std::vector<std::vector<double>>* knot_vectors,
                      const double* control_points,
                      const double* weights);

  virtual constexpr int SplinepyParaDim() const = 0;
  virtual constexpr int SplinepyDim() const = 0;
  virtual std::string SplinepySplineName() const = 0;
  virtual std::string SplinepyWhatAmI() const = 0;
  virtual bool SplinepyHasKnotVectors() const = 0;
  virtual bool SplinepyIsRational() const = 0;
  virtual int SplinepyNumberOfControlPoints() const = 0;
  /// Extract core spline properties. Similar to previous update_p
  virtual void
  SplinepyCurrentProperties(double* degrees,
                            std::vector<std::vector<double>>* knot_vectors,
                            double* control_points,
                            double* weights) const = 0;

  /// Parameter space AABB
  virtual void SplinepyParametricBounds(double* p_bounds) const {
    splinepy::utils::PrintAndThrowError(
        "SplinepyParametricBounds not implemented for",
        SplinepyWhatAmI());
  };

  /// Spline evaluation
  virtual void SplinepyEvaluate(const double* para_coord,
                                double* evaluated) const {
    splinepy::utils::PrintAndThrowError("SplinepyEvaluate not implemented for",
                                        SplinepyWhatAmI());
  };

  /// Spline derivatives
  virtual void SplinepyDerivative(const double* para_coord,
                                  const int* orders,
                                  double* derived) const {
    splinepy::utils::PrintAndThrowError(
        "SplinepyDerivative not implemented for",
        SplinepyWhatAmI());
  };

  /// Spline degree elevation 
  virtual void SplinepyElevateDegree(const int& para_dims) {
    splinepy::utils::PrintAndThrowError(
        "SplinepyElevateDegree not implemented for",
        SplinepyWhatAmI());
  };

  /// Spline degree reduction
  virtual bool SplinepyReduceDegree(const int& para_dims,
                                    const double& tolerance) {
    splinepy::utils::PrintAndThrowError(
        "SplinepyReduceDegree not implemented for",
        SplinepyWhatAmI());
  };

  /// Spline knot insertion.
  virtual void SplinepyInsertKnot(const int& para_dim, const double& knot) {
    splinepy::utils::PrintAndThrowError(
        "SplinepyInsertKnot not implemented for",
        SplinepyWhatAmI());
  };

  /// Spline knot removal.
  virtual bool SplinepyRemoveKnot(const int& para_dim,
                                  const double& knot,
                                  const double& tolerance) {
    splinepy::utils::PrintAndThrowError(
        "SplinepyRemoveKnot not implemented for",
        SplinepyWhatAmI());
  };

  /// Spline multiplication.
  virtual std::shared_ptr<SplinepyBase> SplinepyMultiply(
        const std::shared_ptr<SplinepyBase>& a) const {
     splinepy::utils::PrintAndThrowError(
        "SplinepyMultiply not implemented for",
        SplinepyWhatAmI());
  }

  /// Spline addition.
  virtual std::shared_ptr<SplinepyBase> SplinepyAdd(
    std::shared_ptr<SplinepyBase>& a) const {
     splinepy::utils::PrintAndThrowError(
        "SplinepyAdd not implemented for",
        SplinepyWhatAmI());
  }

  /// Spline composition.
  virtual std::shared_ptr<SplinepyBase> SplinepyCompose(
    std::shared_ptr<SplinepyBase>& inner_function) const {
     splinepy::utils::PrintAndThrowError(
        "SplinepyCompose not implemented for",
        SplinepyWhatAmI());
  }

  /// Check if name matches and throw(=raise) if desired
  virtual bool SplinepySplineNameMatches(const SplinepyBase& a,
                                         const SplinepyBase& b,
                                         const std::string description="",
                                         const bool raise=false) const {
    if (a.SplinepySplineName() != b.SplinepySplineName()) {
      if (raise) {
        splinepy::utils::PrintAndThrowError(
          description,
          "Spline name mismatch -"
          "Spline0:", a.SplinepySplineName(),
          "Spline1:", b.SplinepySplineName()
        );
      }
      return false;
    }
    return true;
  }

  /// Check if para_dim matches and throw(=raise) if desired
  virtual bool SplinepyParaDimMatches(const SplinepyBase& a,
                                      const SplinepyBase& b,
                                      const std::string description="",
                                      const bool raise=false) const {
    if (a.SplinepyParaDim() != b.SplinepyParaDim()) {
      if (raise) {
        splinepy::utils::PrintAndThrowError(
          description,
          "Spline parametric dimension mismatch - "
          "Spline0:", a.SplinepyParaDim(),
          "Spline1:", b.SplinepyParaDim()
        );
      }
      return false;
    }
    return true;
  }

  /// Check if dim matches and throw(=raise) if desired
  virtual bool SplinepyDimMatches(const SplinepyBase& a,
                                  const SplinepyBase& b,
                                  const std::string description="",
                                  const bool raise=false) const {
    if (a.SplinepyDim() != b.SplinepyDim()) {
      if (raise) {
        splinepy::utils::PrintAndThrowError(
          description,
          "Spline parametric dimension mismatch - "
          "Spline0:", a.SplinepyDim(),
          "Spline1:", b.SplinepyDim()
        );
      }
      return false;
    }
    return true;
  }

};

} // namespace splinepy::splines
