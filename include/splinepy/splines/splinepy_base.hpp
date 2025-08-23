/*
MIT License

Copyright (c) 2021 Jaewook Lee

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#pragma once

#include <memory>
#include <vector>

#include "splinepy/utils/arrays.hpp"
#include "splinepy/utils/coordinate_pointers.hpp"

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
  /// Dynamic arrays
  using Array1D_ = splinepy::utils::Array<double, 1>;
  using Array2D_ = splinepy::utils::Array<double, 2>;
  using Array1I_ = splinepy::utils::Array<int, 1>;

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
  /// @brief Extract core spline properties. TODO: should be protected
  /// @param degrees
  /// @param knot_vectors
  /// @param control_points
  /// @param weights
  virtual void
  SplinepyCurrentProperties(int* degrees,
                            std::vector<std::vector<double>>* knot_vectors,
                            double* control_points,
                            double* weights) const = 0;

  virtual std::shared_ptr<ControlPointPointers_> SplinepyControlPointPointers();

  /// @brief Parameter space AABB
  /// @param para_bounds
  virtual void SplinepyParametricBounds(double* para_bounds) const;
  virtual Array2D_ SplinepyParametricBounds() const;

  /// Control mesh resoltuons - number of control points per para dim
  /// @param control_mesh_res
  virtual void SplinepyControlMeshResolutions(int* control_mesh_res) const;
  virtual Array1I_ SplinepyControlMeshResolutions() const;

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
  virtual Array1D_
  SplinepyGrevilleAbscissae(const int& i_para_dim,
                            const double& duplicate_tolerance) const;
  /// @brief Greville abscissae points throughout the whole spline
  /// @param duplicate_tolerance
  /// @return
  virtual Array2D_
  SplinepyGrevilleAbscissae(const double& duplicate_tolerance) const;

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
                                        const bool tight_bounds,
                                        double* para_coord,
                                        double* phys_coord,
                                        double* phys_diff,
                                        double& distance,
                                        double& convergence_norm,
                                        double* first_derivatives,
                                        double* second_derivatives) const;

  /// Spline degree elevation
  virtual void SplinepyElevateDegree(const int& para_dims,
                                     const int multiplicity = 1);

  /// Bezier patch extraction
  virtual std::vector<std::shared_ptr<SplinepyBase>>
  SplinepyExtractBezierPatches() const;

  /// Boundary spline extraction - TODO: const
  virtual std::shared_ptr<SplinepyBase>
  SplinepyExtractBoundary(const int& boundary_id);

  /// Deep copy of current spline
  virtual std::shared_ptr<SplinepyBase> SplinepyDeepCopy() const;
};

} // namespace splinepy::splines
