#ifndef SRC_UTILS_BINOMIALCOEFFICIENTLOOKUPTABLECREATOR_HPP
#define SRC_UTILS_BINOMIALCOEFFICIENTLOOKUPTABLECREATOR_HPP

#include <array>

namespace beziermanipulation::utils {

/*
 * Class to create a look-up table for binomial coefficients
 *
 * @tparam Scalar intrinsic type used to create the lookup table
 * @tparam max_degree maximum depth of the binomial coefficient
 */
template <typename Scalar, unsigned int max_degree>
class BinomialCoefficientLookupTableCreator {
 public:
  static constexpr std::size_t n_values = max_degree * (max_degree + 1) / 2;

  static constexpr std::array<Scalar, n_values> calculate_table() {
    std::array<Scalar, n_values> ret_val{};
    ret_val[0] = static_cast<Scalar>(1);
    for (std::size_t i{0}; i < max_degree; i++) {
      std::size_t offset = i * (i + 1) / 2;
      for (std::size_t j{0}; j <= i; j++) {
        if (i == j || j == 0) {
          ret_val[offset + j] = static_cast<Scalar>(1);
        } else {
          ret_val[offset + j] =
              ret_val[offset - i + j - 1] + ret_val[offset - i + j];
        }
      }
    }
    return ret_val;
  }
};

}

#endif  // SRC_UTILS_BINOMIALCOEFFICIENTLOOKUPTABLECREATOR_HPP
