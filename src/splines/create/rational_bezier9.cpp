#ifdef SPLINEPY_MORE

#include <splinepy/splines/create/create_rational_bezier.hpp>

namespace splinepy::splines::create {

/// dynamic creation of templated Bezier
std::shared_ptr<splinepy::splines::SplinepyBase>
CreateRationalBezier9(const int dim,
                      const int* degrees,
                      const double* control_points,
                      const double* weights) {
  switch (dim) {
  case 1:
    return std::make_shared<RationalBezier<9, 1>>(degrees,
                                                  control_points,
                                                  weights);
  case 2:
    return std::make_shared<RationalBezier<9, 2>>(degrees,
                                                  control_points,
                                                  weights);
  case 3:
    return std::make_shared<RationalBezier<9, 3>>(degrees,
                                                  control_points,
                                                  weights);
  case 4:
    return std::make_shared<RationalBezier<9, 4>>(degrees,
                                                  control_points,
                                                  weights);
  case 5:
    return std::make_shared<RationalBezier<9, 5>>(degrees,
                                                  control_points,
                                                  weights);
  case 6:
    return std::make_shared<RationalBezier<9, 6>>(degrees,
                                                  control_points,
                                                  weights);
  case 7:
    return std::make_shared<RationalBezier<9, 7>>(degrees,
                                                  control_points,
                                                  weights);
  case 8:
    return std::make_shared<RationalBezier<9, 8>>(degrees,
                                                  control_points,
                                                  weights);
  case 9:
    return std::make_shared<RationalBezier<9, 9>>(degrees,
                                                  control_points,
                                                  weights);
  case 10:
    return std::make_shared<RationalBezier<9, 10>>(degrees,
                                                   control_points,
                                                   weights);
  default:
    splinepy::utils::PrintAndThrowError(
        "Something went wrong during CreateBezier. Please help us by writing "
        "an issue about this case at [ github.com/tataratat/splinepy ]");
    break;
  }
  splinepy::utils::PrintAndThrowError(
      "Something went very wrong during CreateBezier. Please help us by "
      "writing "
      "an issue about this case at [ github.com/tataratat/splinepy ]");
  // make compiler happy
  return std::shared_ptr<SplinepyBase>{};
}

} // namespace splinepy::splines::create

#endif
