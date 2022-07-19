/*
MIT License

Copyright (c) 2022 zwar@ilsb.tuwien.ac.at

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

#ifndef SRC_UTILS_BASE64_HPP
#define SRC_UTILS_BASE64_HPP

#include <array>
#include <cassert>
#include <vector>

#include "bezman/src/point.hpp"
#include "bezman/src/utils/type_traits/is_point.hpp"

namespace bezman::utils {

/**
 * @brief Encode for base64 export
 *
 * Python readable base64 export for double values. This allows to store doubles
 * exactly in textformat
 *
 */
template <typename FOO = void*>
class Base64_ {
 private:
  /// Alias for one byte type
  using ByteRepresentation = unsigned char;

  /// Look up table
  static constexpr std::array<char, 64> char_encode_table{
      'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
      'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
      'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
      'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
      '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '+', '/'};

  /**
   * @brief Reverse the encoding table in an array
   *
   * chars are one byte in size
   *
   * @return constexpr std::array<unsigned, 256>
   */
  static constexpr std::array<unsigned, 256> ReverseCharEncodeTable_() {
    static_assert(sizeof(char{}) == 1, "Invalid char size.");
    std::array<unsigned, 256> et_reversed{};
    for (unsigned i{}; i < char_encode_table.size(); i++) {
      et_reversed[static_cast<unsigned>(char_encode_table[i])] = i;
    }
    return et_reversed;
  }

  /// Lookup Table for Decoding B64 string
  static constexpr std::array<unsigned, 256> char_decode_table =
      ReverseCharEncodeTable_();

 public:
  /**
   * @brief Actual encoding routine
   *
   * @tparam BaseType type of individual data entries
   * @param data_vector data to be encoded
   * @return std::string encoded data
   */
  template <
      typename BaseType,
      std::enable_if_t<!type_traits::isPoint_v<BaseType>, void*> = nullptr>
  static std::string Encode(const std::vector<BaseType>& data_vector) {
    const ByteRepresentation* vector_as_bytes =
        reinterpret_cast<const ByteRepresentation*>(&data_vector[0]);

    // Number of bytes for an entry
    constexpr const std::size_t length_of_entry{sizeof(BaseType{})};
    // Minimum number of bytes required
    const std::size_t minimum_n_bytes_required =
        length_of_entry * data_vector.size();
    // Number of bytes must be divisible by three
    const std::size_t additional_padding_bytes =
        (3 - minimum_n_bytes_required % 3) % 3;
    // Required groups of three
    const std::size_t number_of_groups =
        (minimum_n_bytes_required + additional_padding_bytes) / 3;

    // Initialize return value
    std::string encoded_string;
    encoded_string.resize(number_of_groups * 4);

    // Loop over bytes and decode them
    for (std::size_t i_group{}; i_group < number_of_groups; i_group++) {
      const std::size_t buffer_index = i_group * 3;
      std::array<ByteRepresentation, 3> buffer{};
      buffer[0] = buffer_index < minimum_n_bytes_required
                      ? vector_as_bytes[buffer_index + 0]
                      : 0;
      buffer[1] = buffer_index < minimum_n_bytes_required
                      ? vector_as_bytes[buffer_index + 1]
                      : 0;
      buffer[2] = buffer_index < minimum_n_bytes_required
                      ? vector_as_bytes[buffer_index + 2]
                      : 0;

      encoded_string[i_group * 4 + 0] =
          char_encode_table[((buffer[0] & 0xfc) >> 2)];
      encoded_string[i_group * 4 + 1] =
          char_encode_table[((buffer[0] & 0x03) << 4) +
                            ((buffer[1] & 0xf0) >> 4)];
      encoded_string[i_group * 4 + 2] =
          char_encode_table[((buffer[1] & 0x0f) << 2) +
                            ((buffer[2] & 0xc0) >> 6)];
      encoded_string[i_group * 4 + 3] =
          char_encode_table[((buffer[2] & 0x3f) << 0)];
    }

    // Replace trailing invalid data with =
    for (size_t i = 0; i < additional_padding_bytes; ++i) {
      encoded_string[number_of_groups * 4 - i - 1] = '=';
    }

    return encoded_string;
  }

  /**
   * @brief Overload for Point type (handy with splines)
   *
   */
  template <std::size_t dimension, typename Scalar>
  static std::string Encode(
      const std::vector<bezman::Point<dimension, Scalar>>&
          data_vector) {
    std::vector<Scalar> ctps_converted(dimension * data_vector.size());
    for (std::size_t i_point{}; i_point < data_vector.size(); i_point++) {
      for (std::size_t i_dim{}; i_dim < dimension; i_dim++) {
        ctps_converted[i_point * dimension + i_dim] =
            data_vector[i_point][i_dim];
      }
    }
    return Encode(ctps_converted);
  }

  /// Overload for Points
  template <
      typename OutputType,
      std::enable_if_t<type_traits::isPoint_v<OutputType>, void*> = nullptr>
  static std::vector<bezman::Point<OutputType::kSpatialDimension,
                                               typename OutputType::ScalarType>>
  Decode(const std::string& base64string) {
    // Aliases for readability
    using ScalarType = typename OutputType::ScalarType;
    constexpr std::size_t dimension = OutputType::kSpatialDimension;

    // Convert into Scalar Vector
    const std::vector<ScalarType> data_vector =
        Decode<ScalarType>(base64string);

    // Check dimensionality
    const std::size_t n_ctps = data_vector.size() / dimension;
    assert(data_vector.size() % dimension == 0);

    // Init return type
    std::vector<bezman::Point<dimension, ScalarType>> ct_points(
        n_ctps);

    // Transfer data
    for (std::size_t i_point{}; i_point < n_ctps; i_point++) {
      for (std::size_t i_dim{}; i_dim < dimension; i_dim++) {
        ct_points[i_point][i_dim] = data_vector[i_point * dimension + i_dim];
      }
    }
    return ct_points;
  }

  /**
   * @brief Reading a b64 string, transforming it into a vector of a specific
   * type
   *
   * @tparam OutputType target type
   */
  template <
      typename OutputType,
      std::enable_if_t<!type_traits::isPoint_v<OutputType>, void*> = nullptr>
  static std::vector<OutputType> Decode(const std::string& base64string) {
    // Check validity of string
    assert(base64string.size() % 4 == 0);

    // Init return value
    const std::size_t number_of_groups{base64string.size() / 4};
    constexpr const std::size_t length_of_entry{sizeof(OutputType{})};
    const std::size_t number_of_output_values{(number_of_groups * 3) /
                                              length_of_entry};
    std::vector<OutputType> return_value;
    return_value.resize(number_of_output_values);

    // Access as byte stream
    ByteRepresentation* vector_as_bytes =
        reinterpret_cast<ByteRepresentation*>(&return_value[0]);

    // Start the reverse process
    for (std::size_t i_group{}; i_group < number_of_groups; i_group++) {
      const std::size_t buffer_index = i_group * 4;
      std::array<unsigned, 4> buffer{};
      for (unsigned i{}; i < 4; i++) {
        buffer[i] = base64string[buffer_index + i] != '='
                        ? char_decode_table[static_cast<unsigned>(
                              base64string[buffer_index + i])]
                        : 255;
      }

      // Write bytes
      if (buffer[1] != 255) {
        vector_as_bytes[i_group * 3] =
            ((buffer[0] & 0x3f) << 2) + ((buffer[1] & 0x30) >> 4);
      }
      if (buffer[2] != 255) {
        vector_as_bytes[i_group * 3 + 1] =
            ((buffer[1] & 0x0f) << 4) + ((buffer[2] & 0x3c) >> 2);
      }
      if (buffer[3] != 255) {
        vector_as_bytes[i_group * 3 + 2] =
            ((buffer[2] & 0x03) << 6) + ((buffer[3] & 0x3f) >> 0);
      }
    }
    return return_value;
  }
};

/// Non-template implementation can't use compile time functions as constexpr
using Base64 = Base64_<>;

}  // namespace bezman::utils

#endif  // SRC_UTILS_BASE64_HPP
