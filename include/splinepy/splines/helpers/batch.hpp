#include "splinepy/splines/splinepy_base.hpp"
#include "splinepy/utils/array.hpp"

namespace splinepy::splines::helpers {

using SplinepyBase = splinepy::splines::SplinepyBase;

using ArrayI = splinepy::utils::Array<int, 1>;
using ArrayD = splinepy::utils::Array<double, 1>;
using Array2D = splinepy::utils::Array<double, 2>;
using Array2I = splinepy::utils::Array<int, 2>;
using Array3D = splinepy::utils::Array<double, 3>;

Array2D Evaluate(const SplinepyBase& spl,
                 const Array2D& in,
                 double* data_out = nullptr,
                 int n_thread = 1);
Array2D Derivative(const SplinepyBase& spl,
                   const Array2D& in,
                   double* data_out = nullptr,
                   int n_thread = 1);

Array3D Jacobian(const SplinepyBase& spl,
                 const Array2D& in,
                 double* data_out = nullptr,
                 int n_thread = 1);

Array2D JacobianDeterminant(const SplinepyBase& spl,
                            const Array2D& in,
                            double* data_out = nullptr,
                            const int n_thread = 1);

Array2D Basis(const SplinepyBase& spl,
              const Array2D& in,
              double* data_out = nullptr,
              const int n_thread = 1);

Array2I Support(const SplinepyBase& spl,
                const Array2D& in,
                int* data_out = nullptr,
                const int n_thread = 1);

Array2D BasisDerivative(const SplinepyBase& spl,
                        const Array2D& in,
                        const Array2I& order,
                        Array2D& out);
void VerboseProximity(const SplinepyBase&,
                      const Array2D& query,
                      const double tolerance,
                      const int max_iterations,
                      const bool tight_bounds,
                      Array2D& para_coord,
                      Array2D& phys_coord,
                      Array2D& phys_diff,
                      Array2D& distance,
                      s Array2D& convergence_norm,
                      Array2D& first_derivatives,
                      Array3D& second_derivatives);

void BezierCompose(const SplinepyBase&, std::vector<const SplinepyBase*>);
void BezierComposeSensitivities(const SplinepyBase&,
                                std::vector<const SplinepyBase*>);
void BezierCompositionDerivative(const SplinepyBase&,

);
} // namespace splinepy::splines::helpers
