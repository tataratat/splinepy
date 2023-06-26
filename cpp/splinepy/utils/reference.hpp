#pragma once

#include <functional>

namespace splinepy::utils {

/// struct to hold reference.
/// easiler manipulation than reference_wrapper, but better than raw ptr
template<typename T>
struct Reference {
  using Type_ = T;
  /// @brief Reference to the value
  T& value_;
};

} // namespace splinepy::utils
