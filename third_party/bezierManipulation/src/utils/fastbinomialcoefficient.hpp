#ifndef SRC_UTILS_FASTBINOMIALCOEFFICIENT_HPP
#define SRC_UTILS_FASTBINOMIALCOEFFICIENT_HPP

#include <array>

#include "bezierManipulation/src/utils/binomialcoefficientlookuptablecreator.hpp"

#ifndef MAX_BINOMIAL_DEGREE
#define MAX_BINOMIAL_DEGREE 30u
#endif

namespace beziermanipulation::utils {

/*
 * Lookup Class with static binomial coeffient
 *
 * Create a Lookup Table at compile time and provides fast and efficient look up
 * of binomial coefficient. Maximum degree directly correlates with the maximum
 * degree of the bezier spline to be considered.
 */
class FastBinomialCoefficient {
 private:
  /// Sclar Default Type
  using ScalarType = std::size_t;

  /// Size of Lookup table if to large, needs to be removed from stack -> heap
  constexpr static auto n_values =
      BinomialCoefficientLookupTableCreator<ScalarType,
                                            MAX_BINOMIAL_DEGREE>::n_values;

  /// Lookup Table from generator
  constexpr static std::array<ScalarType, n_values> look_up =
      BinomialCoefficientLookupTableCreator<
          ScalarType, MAX_BINOMIAL_DEGREE>::calculate_table();

 public:
  /// Binomial coefficient (careful with asserts in -Ofast mode)
  constexpr static ScalarType choose(const unsigned int n,
                                     const unsigned int i) {
    assert(n < MAX_BINOMIAL_DEGREE);
    return look_up[n * (n + 1) / 2 + i];
  }
};
}  // namespace beziermanipulation::utils

#endif  // SRC_UTILS_FASTBINOMIALCOEFFICIENT_HPP
