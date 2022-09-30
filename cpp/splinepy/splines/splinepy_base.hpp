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
        // return CreateRationalBezier(degrees, control_points, weights);
        return SplinepyCreateRationalBezier();
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
  static std::shared_ptr<SplinepyBase> SplinepyCreateRationalBezier(){};
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
  virtual std::string SplinepyWhatAmI() const = 0;
  virtual int SplinepyNumberOfControlPoints() const = 0;
  /// Extract core spline properties. Similar to previous update_p
  virtual void
  SplinepyCurrentProperties(double* degrees,
                            std::vector<std::vector<double>>* knot_vectors,
                            double* control_points,
                            double* weights) const = 0;

  virtual void SplinepyParametricBounds(double* p_bounds) const {
    splinepy::utils::PrintAndThrowError(
        "SplinepyParametricBounds not implemented for",
        SplinepyWhatAmI());
  };

  virtual void SplinepyEvaluate(const double* para_coord,
                                double* evaluated) const {
    splinepy::utils::PrintAndThrowError("SplinepyEvaluate not implemented for",
                                        SplinepyWhatAmI());
  };

  virtual void SplinepyDerivative(const double* para_coord,
                                  const int* orders,
                                  double* derived) const {
    splinepy::utils::PrintAndThrowError(
        "SplinepyDerivative not implemented for",
        SplinepyWhatAmI());
  };

  virtual void SplinepyElevateDegree(const int& para_dims) {
    splinepy::utils::PrintAndThrowError(
        "SplinepyElevateDegree not implemented for",
        SplinepyWhatAmI());
  };

  virtual bool SplinepyReduceDegree(const int& para_dims,
                                    const double& tolerance) {
    splinepy::utils::PrintAndThrowError(
        "SplinepyReduceDegree not implemented for",
        SplinepyWhatAmI());
  };

  virtual void SplinepyInsertKnot(const int& para_dim, const double& knot) {
    splinepy::utils::PrintAndThrowError(
        "SplinepyInsertKnot not implemented for",
        SplinepyWhatAmI());
  };

  virtual bool SplinepyRemoveKnot(const int& para_dim,
                                  const double& knot,
                                  const double& tolerance) {
    splinepy::utils::PrintAndThrowError(
        "SplinepyRemoveKnot not implemented for",
        SplinepyWhatAmI());
  };
};

} // namespace splinepy::splines
