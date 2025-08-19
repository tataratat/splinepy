#include "splinepy/splines/splinepy_bezier.hpp"

namespace splinepy::splines {

std::shared_ptr<SplinepyBezier> SplinepyBezier::SplinepyMultiply(
    const std::shared_ptr<SplinepyBezier>& a) const {
  splinepy::utils::PrintAndThrowError("SplinepyMultiply not implemented for",
                                      SplinepyWhatAmI());
  return std::shared_ptr<SplinepyBezier>{};
}

std::shared_ptr<SplinepyBezier>
SplinepyBezier::SplinepyAdd(const std::shared_ptr<SplinepyBezier>& a) const {
  splinepy::utils::PrintAndThrowError("SplinepyAdd not implemented for",
                                      SplinepyWhatAmI());
  return std::shared_ptr<SplinepyBezier>{};
}

std::shared_ptr<SplinepyBezier> SplinepyBezier::SplinepyCompose(
    const std::shared_ptr<SplinepyBezier>& inner_function) const {
  splinepy::utils::PrintAndThrowError("SplinepyCompose not implemented for",
                                      SplinepyWhatAmI());
  return std::shared_ptr<SplinepyBezier>{};
}

std::vector<std::shared_ptr<SplinepyBezier>>
SplinepyBezier::SplinepyComposeSensitivities(
    const std::shared_ptr<SplinepyBezier>& inner_function) const {
  splinepy::utils::PrintAndThrowError(
      "SplinepyComposeSensitivities not implemented for",
      SplinepyWhatAmI());
  return std::vector<std::shared_ptr<SplinepyBezier>>{};
}

std::vector<std::shared_ptr<SplinepyBezier>>
SplinepyBezier::SplinepySplit(const int& para_dim,
                              const double& location) const {
  splinepy::utils::PrintAndThrowError("SplinepySplit not implemented for",
                                      SplinepyWhatAmI());
  return {std::shared_ptr<SplinepyBezier>{}};
}

std::shared_ptr<SplinepyBezier>
SplinepyBezier::SplinepyDerivativeSpline(const int* orders) const {
  splinepy::utils::PrintAndThrowError(
      "SplinepyDerivativeSpline is not implemented for",
      SplinepyWhatAmI());
  return std::shared_ptr<SplinepyBezier>{};
}

std::shared_ptr<SplinepyBezier>
SplinepyBezier::SplinepyExtractDim(const int& phys_dim) const {
  splinepy::utils::PrintAndThrowError(
      "SplinepyExtractDim is not implemented for",
      SplinepyWhatAmI());
  return {std::shared_ptr<SplinepyBezier>{}};
}

std::shared_ptr<SplinepyBezier> SplinepyBezier::SplinepyCompositionDerivative(
    const std::shared_ptr<SplinepyBezier>& inner,
    const std::shared_ptr<SplinepyBezier>& inner_derivative) const {
  splinepy::utils::PrintAndThrowError(
      "SplinepyCompositionDerivative is not implemented for",
      SplinepyWhatAmI());
  return {std::shared_ptr<SplinepyBezier>{}};
}
} // namespace splinepy::splines
