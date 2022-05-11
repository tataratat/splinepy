#pragma once

#include <array>
#include <math.h>

#include <Sources/Utilities/named_type.hpp>
#include <Sources/VectorSpaces/vector_space.hpp>
#include <Sources/ParameterSpaces/parameter_space.hpp>
  
using namespace splinelib::sources;  

/* elementwise subtraction */
template<typename T, int dim>
inline void ew_minus(const std::array<T, dim>& arr1,
                     const std::array<T, dim>& arr2,
                     std::array<T, dim>& out) {
  for (int i{0}; i < dim; i++) {
    out[i] = arr1[i] - arr2[i];
  }
}


/* elementwise subtraction overload for SplineLib Coordinate*/
template<typename T, int dim>
inline void ew_minus(
    const vector_spaces::VectorSpace<dim>::Coordinate_ arr1,
    const std::array<T, dim>& arr2,
    std::array<T, dim>& out) {
  for (int i{0}; i < dim; i++) {
    out[i] = arr1[splinelib::Index{i}].Get() - arr2[i];
  }
}


/* elementwise subtraction overload for SplineLib Coordinate*/
template<typename T, int dim>
inline void ew_minus(
    const std::array<T, dim>& arr1,
    const vector_spaces::VectorSpace<dim>::Coordinate_ arr2,
    std::array<T, dim>& out) {
  for (int i{0}; i < dim; i++) {
    out[i] = arr1[i] - arr2[splinelib::Index{i}].Get();
  }
}


/* elementwise subtraction overload for ptr */
template<typename T, int dim>
inline void ew_minus(
    const double* arr1,
    const vector_spaces::VectorSpace<dim>::Coordinate_ arr2,
    std::array<T, dim>& out) {

  for (int i{0}; i < dim; i++) {
    out[i] = arr1[i] - arr2[splinelib::Index{i}].Get();
  }

}


/* elementwise mean */
template<typename T, int dim>
inline void ew_mean(const std::array<T, dim>& arr1,
                    const std::array<T, dim>& arr2,
                    const std::array<T, dim>& out) {
  for (int i{0}; i < dim; i++) {
    out[i] = (arr1[i] + arr2[i]) * .5;
  }

}


template<typename T, int dim>
inline void norm2(std::array<T, dim>& in,
                  double& out) {
  double sqsum = 0.;
  for (int i{0}; i < dim; i++) {
    sqsum += in[i] * in[i];
  }

  out = sqrt(sqsum);
}

/* inplace operation for para coord clipping 
 *
 * Parameters
 * -----------
 * bounds: in
 * para_coord: out
 * clipped: out
 *   -1 -> minimum clip
 *    0 -> no clip
 *    1 -> maximum clip
 */
template<typename T, int para_dim>
inline void clip(
    std::array<std::array<double, para_dim>, 2>& bounds,
    parameter_spaces::ParameterSpace<para_dim>::ParameterCoordinate_& arr1,
    std::array<int, para_dim> clipped) {

  using PC =
    typename parameter_spaces::ParameterSpaces<para_dim>::ParameterCoordinate_;
  using SPC = PC::value_type;

  for (int i{0}; i < para_dim; i++) {
    // check max
    if (arr1[i] > bounds[1][i]) {
      clipped[i] = 1;
      arr1[i] = SPC{bounds[1][i]};
    } else if (arr1[i] < bounds[0][i]) {
      clipped[i] = -1;
      arr1[i] = SPC{bounds[0][i]};
    } else {
      clipped[i] = 0;
    }
  }
}
    
