#include <splinepy/splines/create/create_rational_bezier.hpp>

namespace splinepy::splines::create {

/// dynamic creation of templated Bezier
std::shared_ptr<splinepy::splines::SplinepyBase>
CreateRationalBezier(const int para_dim,
                     const int dim,
                     const int* degrees,
                     const double* control_points,
                     const double* weights) {
  switch (para_dim) {
  case 1:
    return std::make_shared<RationalBezier<1>>(degrees,
                                               control_points,
                                               weights,
                                               dim);
  case 2:
    return std::make_shared<RationalBezier<2>>(degrees,
                                               control_points,
                                               weights,
                                               dim);
  case 3:
    return std::make_shared<RationalBezier<3>>(degrees,
                                               control_points,
                                               weights,
                                               dim);
#ifdef SPLINEPY_MORE
  case 4:
    return std::make_shared<RationalBezier<4>>(degrees,
                                               control_points,
                                               weights,
                                               dim);
  case 5:
    return std::make_shared<RationalBezier<5>>(degrees,
                                               control_points,
                                               weights,
                                               dim);
  case 6:
    return std::make_shared<RationalBezier<6>>(degrees,
                                               control_points,
                                               weights,
                                               dim);
  case 7:
    return std::make_shared<RationalBezier<7>>(degrees,
                                               control_points,
                                               weights,
                                               dim);
  case 8:
    return std::make_shared<RationalBezier<8>>(degrees,
                                               control_points,
                                               weights,
                                               dim);
  case 9:
    return std::make_shared<RationalBezier<9>>(degrees,
                                               control_points,
                                               weights,
                                               dim);
  case 10:
    return std::make_shared<RationalBezier<10>>(degrees,
                                                control_points,
                                                weights,
                                                dim);
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
