#pragma once

#include <memory>
#include <vector>

#include <splinepy/utils/coordinate_pointers.hpp>

namespace bsplinelib::parameter_spaces {
class KnotVector;
}

namespace splinepy::splines {

/// Spline base to enable dynamic use of template splines.
/// Member functions are prepended with "Splinepy".
class SplinepyBase {
public:
  /// Beginning pointers to control points.
  using ControlPointPointers_ = splinepy::utils::ControlPointPointers;
  /// Same type, different alias to emphasize "Weighted"
  using WeightedControlPointPointers_ = splinepy::utils::ControlPointPointers;
  /// Pointers to weights
  using WeightPointers_ = splinepy::utils::WeightPointers;

protected:
  /// each class creates only once and returns shared_ptr second time.
  /// not thread safe for first run.
  std::shared_ptr<ControlPointPointers_> control_point_pointers_ = nullptr;

public:
  /// default ctor
  SplinepyBase() = default;

  /// dtor sets invalid flag to control_point_pointers_ to prevent segfault
  virtual ~SplinepyBase() {
    if (control_point_pointers_) {
      control_point_pointers_->invalid_ = true;
      if (control_point_pointers_->weight_pointers_) {
        control_point_pointers_->weight_pointers_->invalid_ = true;
      }
    }
  };

  /// Dynamically create correct type of spline based on input.
  /// Returned as shared pointer of SplinepyBase
  static std::shared_ptr<SplinepyBase>
  SplinepyCreate(const int para_dim = 0,
                 const int dim = 0,
                 const int* degrees = nullptr,
                 const std::vector<std::vector<double>>* knot_vectors = nullptr,
                 double* control_points = nullptr,
                 double* weights = nullptr);

  /// Dynamic creation of templated bezier
  static std::shared_ptr<SplinepyBase>
  SplinepyCreateBezier(const int para_dim,
                       const int dim,
                       const int* degrees,
                       const double* control_points);

  /// Dynamic creation of templated rational bezier
  static std::shared_ptr<SplinepyBase>
  SplinepyCreateRationalBezier(const int para_dim,
                               const int dim,
                               const int* degrees,
                               const double* control_points,
                               const double* weights);

  /// Dynamic creation of templated bspline
  static std::shared_ptr<SplinepyBase>
  SplinepyCreateBSpline(const int para_dim,
                        const int dim,
                        const int* degrees,
                        const std::vector<std::vector<double>>* knot_vectors,
                        double* control_points);

  /// Dynamic creation of templated nurbs
  static std::shared_ptr<SplinepyBase>
  SplinepyCreateNurbs(const int para_dim,
                      const int dim,
                      const int* degrees,
                      const std::vector<std::vector<double>>* knot_vectors,
                      double* control_points,
                      double* weights);

  /// Check if name matches and throw(=raise) if desired
  static bool SplinepySplineNameMatches(const SplinepyBase& a,
                                        const SplinepyBase& b,
                                        const std::string description = "",
                                        const bool raise = false);

  /// Check if para_dim matches and throw(=raise) if desired
  static bool SplinepyParaDimMatches(const SplinepyBase& a,
                                     const SplinepyBase& b,
                                     const std::string description = "",
                                     const bool raise = false);

  /// Check if dim matches and throw(=raise) if desired
  static bool SplinepyDimMatches(const SplinepyBase& a,
                                 const SplinepyBase& b,
                                 const std::string description = "",
                                 const bool raise = false);

  /// @brief Parametric dimension of spline
  virtual int SplinepyParaDim() const = 0;
  /// @brief Physical dimension of spline
  virtual int SplinepyDim() const = 0;
  /// @brief Returns name of spline
  virtual std::string SplinepySplineName() const = 0;
  /// @brief What am I?
  virtual std::string SplinepyWhatAmI() const = 0;
  /// @brief Returns true iff spline has knot vectors. Bezier splines don’t.
  virtual bool SplinepyHasKnotVectors() const = 0;
  /// @brief Returns true iff spline is rational. NURBS is rational, for
  /// example.
  virtual bool SplinepyIsRational() const = 0;
  /// @brief Get number of control points
  virtual int SplinepyNumberOfControlPoints() const = 0;
  /// @brief Get number of supports
  virtual int SplinepyNumberOfSupports() const = 0;
  /// @brief Returns true iff spline is null-spline
  virtual bool SplinepyIsNull() const { return false; };
  /// @brief Extract core spline properties. Similar to previous update_p
  /// @param degrees
  /// @param knot_vectors
  /// @param control_points
  /// @param weights
  virtual void
  SplinepyCurrentProperties(int* degrees,
                            std::vector<std::vector<double>>* knot_vectors,
                            double* control_points,
                            double* weights) const = 0;

  virtual std::shared_ptr<bsplinelib::parameter_spaces::KnotVector>
  SplinepyKnotVector(const int p_dim);

  virtual std::shared_ptr<ControlPointPointers_> SplinepyControlPointPointers();
  virtual std::shared_ptr<WeightedControlPointPointers_>
  SplinepyWeightedControlPointPointers();
  virtual std::shared_ptr<WeightPointers_> SplinepyWeightPointers();

  /// @brief Parameter space AABB
  /// @param para_bounds
  virtual void SplinepyParametricBounds(double* para_bounds) const;

  /// Control mesh resoltuons - number of control points per para dim
  /// @param control_mesh_res
  virtual void SplinepyControlMeshResolutions(int* control_mesh_res) const;

  /// @brief Calculate Greville abscissae for Spline (required for e.g.
  /// collocation)
  ///
  /// @param[out] greville_abscissae pointer to solution
  /// @param[in] i_para_dim parametric dimension
  /// @param[in] duplicate_tolerance if negative two greville abscissae can be
  ///                                equal, positive tolerance to avoid
  ///                                duplication of greville abscissae. Made to
  ///                                comply with C^(-1) splines. Tolerance
  ///                                represents difference between two greville
  ///                                abscissae for them to be considered equal
  virtual void
  SplinepyGrevilleAbscissae(double* greville_abscissae,
                            const int& i_para_dim,
                            const double& duplicate_tolerance) const;

  /// @brief Evaluate spline
  /// @param[in] para_coord Parametric coordinates
  /// @param[out] evaluated
  virtual void SplinepyEvaluate(const double* para_coord,
                                double* evaluated) const;

  /// @brief Evaluate spline derivatives
  /// @param[in] para_coord Parametric coordinates
  /// @param[in] orders
  /// @param[out] derived
  virtual void SplinepyDerivative(const double* para_coord,
                                  const int* orders,
                                  double* derived) const;

  /// @brief Evaluate jacobians on spline
  /// @param[in] para_coord Parametric coordinates
  /// @param[out] jacobian
  virtual void SplinepyJacobian(const double* para_coord,
                                double* jacobian) const;

  /// @brief Retrieve basis
  /// @param[in] para_coord Parametric coordinates
  /// @param[out] basis
  virtual void SplinepyBasis(const double* para_coord, double* basis) const;

  /// @brief Retrieve basis function derivative
  /// @param[in] para_coord Parametric coordinates
  /// @param[in] order
  /// @param[out] basis
  virtual void SplinepyBasisDerivative(const double* para_coord,
                                       const int* order,
                                       double* basis) const;

  /// Spline Support IDs
  virtual void SplinepySupport(const double* para_coord, int* support) const;

  /// Basis Function values and their support IDs
  virtual void SplinepyBasisAndSupport(const double* para_coord,
                                       double* basis,
                                       int* support) const;

  /// Basis Function Derivative and their support IDs
  virtual void SplinepyBasisDerivativeAndSupport(const double* para_coord,
                                                 const int* orders,
                                                 double* basis,
                                                 int* support) const;

  /// Plants KdTree of sampled spline with given resolution.
  /// KdTree is required for proximity queries.
  virtual void SplinepyPlantNewKdTreeForProximity(const int* resolutions,
                                                  const int& nthreads);

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
                                        double* second_derivatives) const;

  /// Spline degree elevation
  virtual void SplinepyElevateDegree(const int& para_dims);

  /// Spline degree reduction
  virtual bool SplinepyReduceDegree(const int& para_dims,
                                    const double& tolerance);

  /// Spline knot insertion.
  virtual bool SplinepyInsertKnot(const int& para_dim, const double& knot);

  /// Spline knot removal.
  virtual bool SplinepyRemoveKnot(const int& para_dim,
                                  const double& knot,
                                  const double& tolerance);

  /// Spline knot multiplicity per dimension
  virtual std::vector<std::vector<int>> SplinepyKnotMultiplicities() const;

  /// Spline multiplication.
  virtual std::shared_ptr<SplinepyBase>
  SplinepyMultiply(const std::shared_ptr<SplinepyBase>& a) const;

  /// Spline addition.
  virtual std::shared_ptr<SplinepyBase>
  SplinepyAdd(const std::shared_ptr<SplinepyBase>& a) const;

  /// Spline composition.
  virtual std::shared_ptr<SplinepyBase>
  SplinepyCompose(const std::shared_ptr<SplinepyBase>& inner_function) const;

  /// Spline composition sensitivities with respect to the outer spline's
  /// control point positions
  virtual std::vector<std::shared_ptr<SplinepyBase>>
  SplinepyComposeSensitivities(
      const std::shared_ptr<SplinepyBase>& inner_function) const;

  /// Spline Split - single split
  virtual std::vector<std::shared_ptr<SplinepyBase>>
  SplinepySplit(const int& para_dim, const double& location) const;

  /// Derivative spline
  virtual std::shared_ptr<SplinepyBase>
  SplinepyDerivativeSpline(const int* orders) const;

  /// Bezier patch extraction
  virtual std::vector<std::shared_ptr<SplinepyBase>>
  SplinepyExtractBezierPatches() const;

  /// Boundary spline extraction - TODO: const
  virtual std::shared_ptr<SplinepyBase>
  SplinepyExtractBoundary(const int& boundary_id);

  /// Scalar Spline extraction from dim - TODO: const
  virtual std::shared_ptr<SplinepyBase>
  SplinepyExtractDim(const int& phys_dim) const;

  /// Derivative of composition
  virtual std::shared_ptr<SplinepyBase> SplinepyCompositionDerivative(
      const std::shared_ptr<SplinepyBase>& inner,
      const std::shared_ptr<SplinepyBase>& inner_derivative) const;

  /// Deep copy of current spline
  virtual std::shared_ptr<SplinepyBase> SplinepyDeepCopy() const;
};

} // namespace splinepy::splines
