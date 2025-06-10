#include "splinepy/splines/splinepy_base.hpp"
#include "splinepy/utils/array.hpp"

namespace splinepy::splines::helpers {

using SplinepyBase = splinepy::splines::SplinepyBase;

using ArrayI = splinepy::utils::Array<int, 1>;
using ArrayD = splinepy::utils::Array<double, 1>;
using Array2D = splinepy::utils::Array<double, 2>;
using Array2I = splinepy::utils::Array<int, 2>;

void Evaluate(const SplinepyBase& spl,
              const Array2D& query,
              const int n_thread,
              Array2D& out);
void Derivative(const SplinepyBase& spl, const Array2D& query, Array2D& out);

} // namespace splinepy::splines::helpers
