#include <splinepy/splines/create/create_rational_bezier.hpp>

namespace splinepy::splines::create {

/// dynamic creation of templated Bezier
std::shared_ptr<splinepy::splines::SplinepyBase>
CreateRationalBezier1(const int dim,
                      const int* degrees,
                      const double* control_points,
                      const double* weights) {
  switch (dim) {
  case 1:
    return std::make_shared<RationalBezier<1, 1>>(degrees,
                                                  control_points,
                                                  weights);
  case 2:
    return std::make_shared<RationalBezier<1, 2>>(degrees,
                                                  control_points,
                                                  weights);
  case 3:
    return std::make_shared<RationalBezier<1, 3>>(degrees,
                                                  control_points,
                                                  weights);
#ifdef SPLINEPY_MORE
  case 4:
    return std::make_shared<RationalBezier<1, 4>>(degrees,
                                                  control_points,
                                                  weights);
  case 5:
    return std::make_shared<RationalBezier<1, 5>>(degrees,
                                                  control_points,
                                                  weights);
  case 6:
    return std::make_shared<RationalBezier<1, 6>>(degrees,
                                                  control_points,
                                                  weights);
  case 7:
    return std::make_shared<RationalBezier<1, 7>>(degrees,
                                                  control_points,
                                                  weights);
  case 8:
    return std::make_shared<RationalBezier<1, 8>>(degrees,
                                                  control_points,
                                                  weights);
  case 9:
    return std::make_shared<RationalBezier<1, 9>>(degrees,
                                                  control_points,
                                                  weights);
  case 10:
    return std::make_shared<RationalBezier<1, 10>>(degrees,
                                                   control_points,
                                                   weights);
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
