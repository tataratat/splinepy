#include "splinepy/splines/helpers/batched_queries.hpp"
#include "splinepy/utils/nthreads.hpp"
#include "splinepy/utils/print.hpp"

namespace splinepy::splines::helpers {

template<int dim, typename Array>
bool SizeCheck(const Array& arr, const int height, const bool throw_ = true) {
  if (arr.Shape(dim) != height) {
    if (throw_) {
      splinepy::utils::PrintAndThrowError("Expected array size along dim (",
                                          dim,
                                          ") :" height,
                                          "/ current size:",
                                          arr.Shape(dim);
    }
    return false;
  }
  return true;
}

void Evaluate(const SplinepyBase& spline,
              const Array2D& queries,
              Array2D& results,
              int nthreads) {
  SizeCheck<1>(queries, spline.SplinepyParaDim());
  SizeCheck<1>(results, spline.SplinepyDim());
  SizeCheck<0>(results, queries.Shape(0));

  auto evaluate = [&](const int begin, const int end, int) {
    for (int i{begin}; i < end; ++i) {
      spline.SplinepyEvaluate(&queries(i, 0), &results(i, 0));
    }
  };

  splinepy::utils::NThreadExecution(evaluate, queries.Shape(0), nthreads);
}

void Derivative(const SplinepyBase& spline,
                const Array2D& queries,
                const Array2I& orders,
                Array2D& results,
                int nthreads) {
  const int para_dim = spline.SplinepyParaDim();
  const int dim = spline.SplinepyDim();

  SizeCheck<1>(queries, para_dim);
  SizeCheck<1>(orders, para_dim);
  SizeCheck<1>(results, dim);
  SizeCheck<0>(results, queries.Shape(0) * orders.Shape(0));

  auto derivative = [&](const int begin, const int end, int) {
    const int outer_stride = orders.Shape(0);
    for (int i{begin}; i < end; ++i) {
      for (int j{}; j < orders.Shape(0); ++j) {
        spline.SplinepyDerivative(&queries(i, 0),
                                  &orders(j, 0),
                                  &results(i * outer_stride + j, 0));
      }
    }
  };

  splinepy::utils::NThreadExecution(derivative, queries.Shape(0), nthreads);
}
void UniformSample(const SplinepyBase& spline,
                   const Array1I& resolutions,
                   Array2D& results,
                   int nthreads);
void UniformSampleDerivative(const SplinepyBase& spline,
                             const Array1I& resolutions,
                             const Array2I& orders,
                             Array2D& results,
                             int nthreads);
void Support(const SplinepyBase& spline,
             const Array2D& queries,
             Array2I& support,
             int nthreads);
void Basis(const SplinepyBase& spline,
           const Array2D& queries,
           Array2D& basis,
           int nthreads);
void BasisDerivative(const SplinepyBase& spline,
                     const Array2D& queries,
                     Array2D& basis_derivative,
                     int nthreads);
void Proximities(const SplinepyBase& spline,
                 const Array2D& queries,
                 Array2D& para_coords,
                 Array2D& phys_coords,
                 Array2D& phys_diff,
                 Array1D& distance,
                 Array1D& convergence_norm,
                 Array3D& first_derivatives,
                 Array4D& second_derivatives);
} // namespace splinepy::splines::helpers
