#include <splinepy/splines/create/create_bezier.hpp>

namespace splinepy::splines::create {

/// dynamic creation of templated Bezier
std::shared_ptr<splinepy::splines::SplinepyBase>
CreateBezier(const int para_dim,
             const int dim,
             const int* degrees,
             const double* control_points) {
  switch (para_dim) {
  case 1:
    return std::make_shared<Bezier<1>>(degrees, control_points, dim);
  case 2:
    return std::make_shared<Bezier<2>>(degrees, control_points, dim);
  case 3:
    return std::make_shared<Bezier<3>>(degrees, control_points, dim);
#ifdef SPLINEPY_MORE
  case 4:
    return std::make_shared<Bezier<4>>(degrees, control_points, dim);
  case 5:
    return std::make_shared<Bezier<5>>(degrees, control_points, dim);
  case 6:
    return std::make_shared<Bezier<6>>(degrees, control_points, dim);
  case 7:
    return std::make_shared<Bezier<7>>(degrees, control_points, dim);
  case 8:
    return std::make_shared<Bezier<8>>(degrees, control_points, dim);
  case 9:
    return std::make_shared<Bezier<9>>(degrees, control_points, dim);
  case 10:
    return std::make_shared<Bezier<10>>(degrees, control_points, dim);
#endif
  default:
    splinepy::utils::PrintAndThrowError(
        "Something went wrong during CreateBezier. Please help us by writing "
        "an issue about this case at [ github.com/tataratat/splinepy ]");
  }
  splinepy::utils::PrintAndThrowError(
      "Something went very wrong during CreateBezier. Please help us by "
      "writing "
      "an issue about this case at [ github.com/tataratat/splinepy ]");
  // make compiler happy
  return std::shared_ptr<SplinepyBase>{};
}

} // namespace splinepy::splines::create
