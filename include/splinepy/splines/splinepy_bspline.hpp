#pragma once

#include "splinepy/splines/splinepy_base.hpp"

namespace bsplinelib::parameter_spaces {
class KnotVector;
class ParameterSpaceBase;
} // namespace bsplinelib::parameter_spaces

namespace splinepy::splines {

class SplinepyBSpline : public splinepy::splines::SplinepyBase {
public:
  virtual std::shared_ptr<bsplinelib::parameter_spaces::ParameterSpaceBase>
  SplinepyParameterSpace();
  virtual std::shared_ptr<bsplinelib::parameter_spaces::KnotVector>
  SplinepyKnotVector(const int p_dim);

  /// Spline degree reduction
  virtual bool SplinepyReduceDegree(const int& para_dims,
                                    const double& tolerance);
  /// Spline knot insertion.
  virtual int SplinepyInsertKnot(const int& para_dim,
                                 const double& knot,
                                 const int multiplicity = 1);

  /// Spline knot removal.
  virtual bool SplinepyRemoveKnot(const int& para_dim,
                                  const double& knot,
                                  const double& tolerance);

  /// Spline knot multiplicity per dimension
  virtual std::vector<std::vector<int>> SplinepyKnotMultiplicities() const;
};
} // namespace splinepy::splines
