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
  virtual int SplinepyNumberOfSupports() const = 0;
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

  /// Basis Function values and their support IDs
  virtual void SplinepyBasisAndSupport(const double* para_coord,
                                       double* basis,
                                       int* support) const {
    splinepy::utils::PrintAndThrowError(
        "SplinepyBasisAndSupport not implemented for",
        SplinepyWhatAmI());
  }

  virtual void SplinepyPlantNewKdTreeForProximity(const int* resolutions,
                                                  const int& nthreads) {
    splinepy::utils::PrintAndThrowError(
        "SplinepyPlantNewKdTreeForProximity not implemented for",
        SplinepyWhatAmI());
  }

  /// Verbose proximity query - make sure to plant a kdtree first.
  virtual void SplinepyVerboseProximity(const double* query,
                                        const double& tolerance,
                                        const int& max_iterations,
                                        const bool aggressive_bounds,
                                        double* para_coord,
                                        double* phys_coord,
                                        double* phys_diff,
                                        double& distance,
                                        double& convergence_norm,
                                        double* first_derivatives,
                                        double* second_derivatives) const {
    splinepy::utils::PrintAndThrowError(
        "SplinepyVerboseProximity not implemented for",
        SplinepyWhatAmI());
  }

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
    return false;
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
    return false;
  };

  /// Spline multiplication.
  virtual std::shared_ptr<SplinepyBase>
  SplinepyMultiply(const std::shared_ptr<SplinepyBase>& a) const {
    splinepy::utils::PrintAndThrowError("SplinepyMultiply not implemented for",
                                        SplinepyWhatAmI());
    return std::shared_ptr<SplinepyBase>{};
  }

  /// Spline addition.
  virtual std::shared_ptr<SplinepyBase>
  SplinepyAdd(const std::shared_ptr<SplinepyBase>& a) const {
    splinepy::utils::PrintAndThrowError("SplinepyAdd not implemented for",
                                        SplinepyWhatAmI());
    return std::shared_ptr<SplinepyBase>{};
  }

  /// Spline composition.
  virtual std::shared_ptr<SplinepyBase>
  SplinepyCompose(const std::shared_ptr<SplinepyBase>& inner_function) const {
    splinepy::utils::PrintAndThrowError("SplinepyCompose not implemented for",
                                        SplinepyWhatAmI());
    return std::shared_ptr<SplinepyBase>{};
  }

  /// Spline Split - single split
  virtual std::vector<std::shared_ptr<SplinepyBase>>
  SplinepySplit(const int& para_dim, const double& location) const {
    splinepy::utils::PrintAndThrowError("SplinepySplit not implemented for",
                                        SplinepyWhatAmI());
    return {std::shared_ptr<SplinepyBase>{}};
  }

  /// Derivative spline
  virtual std::shared_ptr<SplinepyBase>
  SplinepyDerivativeSpline(const int* orders) const {
    splinepy::utils::PrintAndThrowError(
        "SplinepyDerivativeSpline is not implemented for",
        SplinepyWhatAmI());
    return std::shared_ptr<SplinepyBase>{};
  };

  /// Bezier patch extraction
  virtual std::vector<std::shared_ptr<SplinepyBase>>
  SplinepyExtractBezierPatches() const {
    splinepy::utils::PrintAndThrowError(
        "SplinepyBezierPatchExtraction is not implemented for",
        SplinepyWhatAmI());
    return {std::shared_ptr<SplinepyBase>{}};
  }

  /// Check if name matches and throw(=raise) if desired
  static bool SplinepySplineNameMatches(const SplinepyBase& a,
                                        const SplinepyBase& b,
                                        const std::string description = "",
                                        const bool raise = false) {
    if (a.SplinepySplineName() != b.SplinepySplineName()) {
      if (raise) {
        splinepy::utils::PrintAndThrowError(description,
                                            "Spline name mismatch -"
                                            "Spline0:",
                                            "/",
                                            a.SplinepySplineName(),
                                            "Spline1:",
                                            b.SplinepySplineName());
      }
      return false;
    }
    return true;
  }

  /// Check if para_dim matches and throw(=raise) if desired
  static bool SplinepyParaDimMatches(const SplinepyBase& a,
                                     const SplinepyBase& b,
                                     const std::string description = "",
                                     const bool raise = false) {
    if (a.SplinepyParaDim() != b.SplinepyParaDim()) {
      if (raise) {
        splinepy::utils::PrintAndThrowError(
            description,
            "Spline parametric dimension mismatch - "
            "Spline0:",
            "/",
            a.SplinepyParaDim(),
            "Spline1:",
            b.SplinepyParaDim());
      }
      return false;
    }
    return true;
  }

  /// Check if dim matches and throw(=raise) if desired
  static bool SplinepyDimMatches(const SplinepyBase& a,
                                 const SplinepyBase& b,
                                 const std::string description = "",
                                 const bool raise = false) {
    if (a.SplinepyDim() != b.SplinepyDim()) {
      if (raise) {
        splinepy::utils::PrintAndThrowError(
            description,
            "Spline parametric dimension mismatch - "
            "Spline0:",
            "/",
            a.SplinepyDim(),
            "Spline1:",
            b.SplinepyDim());
      }
      return false;
    }
    return true;
  }
};

} // namespace splinepy::splines
