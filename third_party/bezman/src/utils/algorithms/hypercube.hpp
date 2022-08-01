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

#ifndef UTILS_ALGORITHMS_HYPERCUBE_HPP
#define UTILS_ALGORITHMS_HYPERCUBE_HPP

#include <array>

#include "bezman/src/utils/algorithms/int_power.hpp"

namespace bezman::utils::algorithms {

/**
 * @brief Provides functions and auxiliary values for hypercubes
 *
 * This can be used to retrieve indices for values in the parametric space
 *
 * Element Numberings:
 * ------------------------------------------------
 * 2D:
 *
 *  3 --(2)>- 2
 *  ^         ^
 * (3)   0   (1)
 *  |         |
 *  0 --(0)>- 1
 *
 * directions:
 * x : 0 -> 1
 * y : 0 -> 3
 * ------------------------------------------------
 *
 * 3D:
 *
 *             7 --(2)>- 6
 *             ^         ^
 *            (11)  2   (10)
 *             |         |
 *   7 -<(11)- 3 --(2)>- 2 --(10)> 6 -<(6)-- 7
 *   ^         ^         ^         ^         ^
 *  (7)   3   (3)   0   (1)   1   (5)   5   (7)
 *   |         |         |         |         |
 *   4 -<(8)-- 0 --(0)>- 1 --(9)>- 5 -<(4)-- 4
 *             |         |
 *            (8)   4   (9)
 *             v         v
 *             4 --(0)>- 5
 *
 * directions:
 * x : 0 -> 1
 * y : 0 -> 3
 * z : 0 -> 4
 * ------------------------------------------------
 */
template <std::size_t dimension>
class HyperCube {
  // Current implementation limited to 2 and 3 dimensions
  static_assert((dimension == 3 || dimension == 2),
                "High-Dimensional and Line Patches not supported");

 public:
  /**
   * @brief Get the Opposite Faces
   *
   * See Pictures for numbering system in class header
   */
  static constexpr std::array<std::size_t, dimension * 2> GetOppositeFaces() {
    static_assert((dimension == 3 || dimension == 2),
                  "High-Dimensional and Line Patches not supported");
    if constexpr (dimension == 2) {
      return std::array<std::size_t, dimension * 2>{2, 3, 0, 1};
    } else if constexpr (dimension == 3) {
      // Must be 3D
      return std::array<std::size_t, dimension * 2>{5, 3, 4, 1, 2, 0};
    } else {
      // Should never happen
      return std::array<std::size_t, dimension * 2>{};
    }
  }

  /**
   * @brief Control point indices to corner-vertices
   *
   * The function returns the associated indices in the controlpoint vector with
   * respect to the spline degrees. The enumeration is based on the mfem element
   * definitions
   *
   * See Pictures for numbering system in class header
   */
  static std::array<std::size_t,
                    bezman::utils::algorithms::IntPower(
                        static_cast<std::size_t>(2), dimension)>
  VertexIdForDegrees(const std::array<std::size_t, dimension> &degrees) {
    // Alias for Readability
    using ReturnType =
        std::array<std::size_t, bezman::utils::algorithms::IntPower(
                                    static_cast<std::size_t>(2), dimension)>;
    static_assert((dimension == 3 || dimension == 2),
                  "High-Dimensional and Line Patches not supported");
    if constexpr (dimension == 2) {
      return ReturnType{
          0,                                        // 0
          degrees[0],                               // 1
          (degrees[0] + 1) * (degrees[1] + 1) - 1,  // 2
          (degrees[0] + 1) * degrees[1]             // 3
      };
    } else {
      return ReturnType{
          // bottom: 0 - 1 - 2 - 3
          0,                                        // 0
          degrees[0],                               // 1
          (degrees[0] + 1) * (degrees[1] + 1) - 1,  // 2
          (degrees[0] + 1) * degrees[1],            // 3
          // top: 4 - 5 - 6 - 7
          degrees[2] * (degrees[0] + 1) * (degrees[1] + 1),               // 4
          degrees[2] * (degrees[0] + 1) * (degrees[1] + 1) + degrees[0],  // 5
          (degrees[2] + 1) * (degrees[0] + 1) * (degrees[1] + 1) - 1,     // 6
          degrees[2] * (degrees[0] + 1) * (degrees[1] + 1) +
              (degrees[0] + 1) * degrees[1]  // 7
      };
    }
  }

  /**
   * @brief Local Element Vertex IDs to the associated subelement
   *
   * Returns local element vertices to the associated subelement, i.e. lines for
   * squares, squares for hexahedra, etc...
   *
   * Here it is important, that opposite faces, have the same numbering. I.e.
   * face[i] lies on opposite_face[i] if the elements are neighbors
   *
   * See Pictures for numbering system in class header
   */
  static constexpr std::array<
      std::array<std::size_t, bezman::utils::algorithms::IntPower(
                                  static_cast<std::size_t>(2), dimension - 1)>,
      dimension * 2>
  SubElementVerticesToFace() {
    using ReturnType =
        std::array<std::array<std::size_t,
                              bezman::utils::algorithms::IntPower(
                                  static_cast<std::size_t>(2), dimension - 1)>,
                   2 * dimension>;

    static_assert(dimension == 2 || dimension == 3, "Not Implemented");
    if constexpr (dimension == 2) {
      return ReturnType{0, 1, 1, 2, 3, 2, 0, 3};
    } else if constexpr (dimension == 3) {
      return ReturnType{
          0, 1, 2, 3,  // 0
          1, 2, 6, 5,  // 1
          3, 2, 6, 7,  // 2
          0, 3, 7, 4,  // 3
          0, 1, 5, 4,  // 4
          4, 5, 6, 7   // 5
      };
    }
  }

  /**
   * @brief Local Element Vertex indices associated to the edges
   *
   * See Pictures for numbering system in class header
   * TODO changed this code because of an  compile error
   */
  static constexpr std::array<
      std::array<std::size_t, static_cast<std::size_t>(2)>,
      dimension == 2 ? static_cast<std::size_t>(4)
                     : static_cast<std::size_t>(12)>
  EdgeVertexIndices() {
    using ReturnType =
        std::array<std::array<std::size_t, 2>,
                   dimension == 2 ? static_cast<std::size_t>(4)
                                  : static_cast<std::size_t>(12)>;

    static_assert(dimension == 2 || dimension == 3, "Not Implemented");
    if constexpr (dimension == 2) {
      return ReturnType{// x- direction
                        0, 1, 3, 2,
                        // y-direction
                        1, 2, 0, 3};
    } else if constexpr (dimension == 3) {
      return ReturnType{// Edges in x-direction
                        0, 1, 3, 2, 4, 5, 7, 6,
                        // Edges in y-direction
                        1, 2, 0, 3, 5, 6, 4, 7,
                        // Edgegs in z-direction
                        0, 4, 1, 5, 2, 6, 3, 7};
    }
  }

  /**
   * @brief Returns the indices of the faces, that are normal to the i-th
   * parametric dimension
   *
   * See Pictures for numbering system in class header
   *
   * @return A two dimensional array containing the indices of the faces which
   * are normal to the i-th parametric dimension
   */
  static constexpr std::array<std::array<std::size_t, (2 * dimension) - 2>,
                              dimension>
  GetNormalFaceIndicesToParametricDimension() {
    static_assert(dimension == 2 || dimension == 3, "Not Implemented");
    using ReturnType =
        std::array<std::array<std::size_t, (2 * dimension) - 2>, dimension>;
    if constexpr (dimension == 2) {
      return ReturnType{// X-direction
                        0, 2,
                        // Y-direction
                        1, 3};
    } else if constexpr (dimension == 3) {
      return ReturnType{// x-direction
                        2, 4, 0, 5,
                        // y-direction
                        1, 3, 0, 5,
                        // z-direction
                        1, 3, 2, 4};
    }
  }
};
}  // namespace bezman::utils::algorithms
#endif  // UTILS_ALGORITHMS_HYPERCUBE_HPP
