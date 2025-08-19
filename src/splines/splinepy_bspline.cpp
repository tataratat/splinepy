#include "splinepy/splines/splinepy_bspline.hpp"

namespace splinepy::splines {

std::shared_ptr<bsplinelib::parameter_spaces::ParameterSpaceBase>
SplinepyBSpline::SplinepyParameterSpace() {
  splinepy::utils::PrintAndThrowError(
      "SplinepyParameterSpace not implemented for",
      SplinepyWhatAmI());
  return nullptr;
}

std::shared_ptr<bsplinelib::parameter_spaces::KnotVector>
SplinepyBSpline::SplinepyKnotVector(const int p_dim) {
  splinepy::utils::PrintAndThrowError("SplinepyKnotVector not implemented for",
                                      SplinepyWhatAmI());
  return nullptr;
}

bool SplinepyBSpline::SplinepyReduceDegree(const int& para_dims,
                                           const double& tolerance) {
  splinepy::utils::PrintAndThrowError(
      "SplinepyReduceDegree not implemented for",
      SplinepyWhatAmI());
  return false;
}

int SplinepyBSpline::SplinepyInsertKnot(const int& para_dim,
                                        const double& knot,
                                        const int multiplicity) {
  splinepy::utils::PrintAndThrowError("SplinepyInsertKnot not implemented for",
                                      SplinepyWhatAmI());
  return -1;
}

bool SplinepyBSpline::SplinepyRemoveKnot(const int& para_dim,
                                         const double& knot,
                                         const double& tolerance) {
  splinepy::utils::PrintAndThrowError("SplinepyRemoveKnot not implemented for",
                                      SplinepyWhatAmI());
  return false;
}

std::vector<std::vector<int>>
SplinepyBSpline::SplinepyKnotMultiplicities() const {
  splinepy::utils::PrintAndThrowError(
      "SplinepyKnotMultiplicities not implemented for",
      SplinepyWhatAmI());
  return std::vector<std::vector<int>>{};
};

} // namespace splinepy::splines
