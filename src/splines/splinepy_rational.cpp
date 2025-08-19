#include "splinepy/splines/splinepy_rational.hpp"

namespace splinepy::splines {

std::shared_ptr<typename SplinepyRational::WeightedControlPointPointers_>
SplinepyRational::SplinepyWeightedControlPointPointers() {
  splinepy::utils::PrintAndThrowError(
      "SplinepyWeightedControlPointPointers not implemented for",
      SplinepyWhatAmI());
  return nullptr;
}

std::shared_ptr<typename SplinepyRational::WeightPointers_>
SplinepyRational::SplinepyWeightPointers() {
  splinepy::utils::PrintAndThrowError(
      "SplinepyWeightPointers not implemented for",
      SplinepyWhatAmI());
  return nullptr;
}

} // namespace splinepy::splines
