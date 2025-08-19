#pragma once

#include "splinepy/splines/splinepy_base.hpp"

namespace splinepy::splines {

class SplinepyBezier : public splinepy::splines::SplinepyBase {
public:
  using Base_ = splinepy::splines::SplinepyBase;

  /// Spline multiplication.
  virtual std::shared_ptr<SplinepyBezier>
  SplinepyMultiply(const std::shared_ptr<SplinepyBezier>& a) const;

  /// Spline addition.
  virtual std::shared_ptr<SplinepyBezier>
  SplinepyAdd(const std::shared_ptr<SplinepyBezier>& a) const;

  /// Spline composition.
  virtual std::shared_ptr<SplinepyBezier>
  SplinepyCompose(const std::shared_ptr<SplinepyBezier>& inner_function) const;

  /// Spline composition sensitivities with respect to the outer spline's
  /// control point positions
  virtual std::vector<std::shared_ptr<SplinepyBezier>>
  SplinepyComposeSensitivities(
      const std::shared_ptr<SplinepyBezier>& inner_function) const;

  /// Spline Split - single split
  virtual std::vector<std::shared_ptr<SplinepyBezier>>
  SplinepySplit(const int& para_dim, const double& location) const;

  /// Derivative spline
  virtual std::shared_ptr<SplinepyBezier>
  SplinepyDerivativeSpline(const int* orders) const;

  /// Scalar Spline extraction from dim
  virtual std::shared_ptr<SplinepyBezier>
  SplinepyExtractDim(const int& phys_dim) const;

  /// Derivative of composition
  virtual std::shared_ptr<SplinepyBezier> SplinepyCompositionDerivative(
      const std::shared_ptr<SplinepyBezier>& inner,
      const std::shared_ptr<SplinepyBezier>& inner_derivative) const;
};
} // namespace splinepy::splines
