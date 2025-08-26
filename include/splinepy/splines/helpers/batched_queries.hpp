/// Batched execution of queries. Optionally using multithreading
/// Note that without binding or

#include "splinepy/splines/splinepy_base.hpp"
#include "splinepy/utils/arrays.hpp"

namespace splinepy::splines::helpers {

using SplinepyBase = splinepy::splines::SplinepyBase;
using Array1I = splinepy::utils::Array1I;
using Array2I = splinepy::utils::Array2I;
using Array1D = splinepy::utils::Array1D;
using Array2D = splinepy::utils::Array2D;
using Array3D = splinepy::utils::Array3D;
using Array4D = splinepy::utils::Array4D;

void Evaluate(const SplinepyBase& spline,
              const Array2D& queries,
              Array2D& results,
              int nthreads);
void Derivative(const SplinepyBase& spline,
                const Array2D& queries,
                const Array2I& orders,
                Array2D& results,
                int nthreads);
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
