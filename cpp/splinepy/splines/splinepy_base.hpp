#pragma once

#include <memory>
#include <vector>

#include <splinepy/utils/reference.hpp>

namespace splinepy::splines {

/// Spline base to enable dynamic use of template splines.
/// Member functions are prepended with "Splinepy".
class SplinepyBase {
public:
  using CoordinateReferences_ = std::vector<splinepy::utils::Reference<double>>;

  /// default ctor
  SplinepyBase() = default;
  ///
  virtual ~SplinepyBase(){};

  /// Dynamically create correct type of spline based on input.
  /// Returned as shared pointer of SplinepyBase
  static std::shared_ptr<SplinepyBase>
  SplinepyCreate(const int para_dim = 0,
                 const int dim = 0,
                 const int* degrees = nullptr,
                 const std::vector<std::vector<double>>* knot_vectors = nullptr,
                 const double* control_points = nullptr,
                 const double* weights = nullptr);

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
                        const double* control_points);

  /// Dynamic creation of templated nurbs
  static std::shared_ptr<SplinepyBase>
  SplinepyCreateNurbs(const int para_dim,
                      const int dim,
                      const int* degrees,
                      const std::vector<std::vector<double>>* knot_vectors,
                      const double* control_points,
                      const double* weights);

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

  virtual int SplinepyParaDim() const = 0;
  virtual int SplinepyDim() const = 0;
  virtual std::string SplinepySplineName() const = 0;
  virtual std::string SplinepyWhatAmI() const = 0;
  virtual bool SplinepyHasKnotVectors() const = 0;
  virtual bool SplinepyIsRational() const = 0;
  virtual int SplinepyNumberOfControlPoints() const = 0;
  virtual int SplinepyNumberOfSupports() const = 0;
  /// Extract core spline properties. Similar to previous update_p
  virtual void
  SplinepyCurrentProperties(int* degrees,
                            std::vector<std::vector<double>>* knot_vectors,
                            double* control_points,
                            double* weights) const = 0;

  virtual std::shared_ptr<CoordinateReferences_> SplinepyCoordinateReferences();

  /// Parameter space AABB
  virtual void SplinepyParametricBounds(double* p_bounds) const;

  /// Control mesh resoltuons - number of control points per para dim
  virtual void SplinepyControlMeshResolutions(int* control_mesh_res) const;

  /**
   * @brief Calculate greville abscissae for Spline (required for e.g.
   * collocation)
   *
   * @param greville_abscissae[out] pointer to solution
   * @param i_para_dim[in] parametric dimension
   */
  virtual void SplinepyGrevilleAbscissae(double* greville_abscissae,
                                         const int& i_para_dim) const;

  /// Spline evaluation
  virtual void SplinepyEvaluate(const double* para_coord,
                                double* evaluated) const;

  /// Spline derivatives
  virtual void SplinepyDerivative(const double* para_coord,
                                  const int* orders,
                                  double* derived) const;

  /// Spline evaluation
  virtual void SplinepyJacobian(const double* para_coord,
                                double* jacobian) const;

  /// Basis Function values
  virtual void SplinepyBasis(const double* para_coord, double* basis) const;

  /// Basis Function derivative values
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
  virtual void SplinepyInsertKnot(const int& para_dim, const double& knot);

  /// Spline knot removal.
  virtual bool SplinepyRemoveKnot(const int& para_dim,
                                  const double& knot,
                                  const double& tolerance);

  /// Spline multiplication.
  virtual std::shared_ptr<SplinepyBase>
  SplinepyMultiply(const std::shared_ptr<SplinepyBase>& a) const;

  /// Spline addition.
  virtual std::shared_ptr<SplinepyBase>
  SplinepyAdd(const std::shared_ptr<SplinepyBase>& a) const;

  /// Spline composition.
  virtual std::shared_ptr<SplinepyBase>
  SplinepyCompose(const std::shared_ptr<SplinepyBase>& inner_function) const;

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
};

} // namespace splinepy::splines
