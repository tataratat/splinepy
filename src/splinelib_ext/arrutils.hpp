#pragma once

#include <array>
#include <math.h>

#include <Sources/Utilities/named_type.hpp>
#include <Sources/VectorSpace/vector_space.hpp>

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
    const
        splinelib::sources::vector_spaces::VectorSpace<dim>::Coordinate_ arr1,
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
    const
        splinelib::sources::vector_spaces::VectorSpace<dim>::Coordinate_ arr2,
    std::array<T, dim>& out) {
  for (int i{0}; i < dim; i++) {
    out[i] = arr1[i] - arr2[splinelib::Index{i}].Get();
  }
}


/* elementwise subtraction overload for ptr */
template<typename T, int dim>
inline void ew_minus(
    const double* arr1,
    const
        splinelib::sources::vector_spaces::VectorSpace<dim>::Coordinate_ arr2,
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
