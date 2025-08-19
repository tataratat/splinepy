#pragma once

#include "splinepy/splines/splinepy_base.hpp"

namespace splinepy::splines {

/// @brief Rational splines interface.
class SplinepyRational : public SplinepyBase {
public:
  using Base_ = splinepy::splines::SplinepyBase;
  using WeightedControlPointPointers_ = Base_::WeightedControlPointPointers_;

  virtual std::shared_ptr<WeightedControlPointPointers_>
  SplinepyWeightedControlPointPointers();
  virtual std::shared_ptr<WeightPointers_> SplinepyWeightPointers();
};
} // namespace splinepy::splines
