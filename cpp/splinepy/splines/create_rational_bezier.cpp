#include <splinepy/splines/rational_bezier.hpp>
#include <splinepy/splines/splinepy_base.hpp>

namespace splinepy::splines {

/// dynamic creation of templated RationalBezier
std::shared_ptr<SplinepyBase>
splinepy::splines::SplinepyBase::SplinepyCreateRationalBezier(
    const int para_dim,
    const int dim,
    const int* degrees,
    const double* control_points,
    const double* weights) {
  switch (para_dim) {
  case 1:
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
          "Something went wrong during CreateRationalBezier. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
  case 2:
    switch (dim) {
    case 1:
      return std::make_shared<RationalBezier<2, 1>>(degrees,
                                                    control_points,
                                                    weights);
    case 2:
      return std::make_shared<RationalBezier<2, 2>>(degrees,
                                                    control_points,
                                                    weights);
    case 3:
      return std::make_shared<RationalBezier<2, 3>>(degrees,
                                                    control_points,
                                                    weights);
#ifdef SPLINEPY_MORE
    case 4:
      return std::make_shared<RationalBezier<2, 4>>(degrees,
                                                    control_points,
                                                    weights);
    case 5:
      return std::make_shared<RationalBezier<2, 5>>(degrees,
                                                    control_points,
                                                    weights);
    case 6:
      return std::make_shared<RationalBezier<2, 6>>(degrees,
                                                    control_points,
                                                    weights);
    case 7:
      return std::make_shared<RationalBezier<2, 7>>(degrees,
                                                    control_points,
                                                    weights);
    case 8:
      return std::make_shared<RationalBezier<2, 8>>(degrees,
                                                    control_points,
                                                    weights);
    case 9:
      return std::make_shared<RationalBezier<2, 9>>(degrees,
                                                    control_points,
                                                    weights);
    case 10:
      return std::make_shared<RationalBezier<2, 10>>(degrees,
                                                     control_points,
                                                     weights);
#endif
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateRationalBezier. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
  case 3:
    switch (dim) {
    case 1:
      return std::make_shared<RationalBezier<3, 1>>(degrees,
                                                    control_points,
                                                    weights);
    case 2:
      return std::make_shared<RationalBezier<3, 2>>(degrees,
                                                    control_points,
                                                    weights);
    case 3:
      return std::make_shared<RationalBezier<3, 3>>(degrees,
                                                    control_points,
                                                    weights);
#ifdef SPLINEPY_MORE
    case 4:
      return std::make_shared<RationalBezier<3, 4>>(degrees,
                                                    control_points,
                                                    weights);
    case 5:
      return std::make_shared<RationalBezier<3, 5>>(degrees,
                                                    control_points,
                                                    weights);
    case 6:
      return std::make_shared<RationalBezier<3, 6>>(degrees,
                                                    control_points,
                                                    weights);
    case 7:
      return std::make_shared<RationalBezier<3, 7>>(degrees,
                                                    control_points,
                                                    weights);
    case 8:
      return std::make_shared<RationalBezier<3, 8>>(degrees,
                                                    control_points,
                                                    weights);
    case 9:
      return std::make_shared<RationalBezier<3, 9>>(degrees,
                                                    control_points,
                                                    weights);
    case 10:
      return std::make_shared<RationalBezier<3, 10>>(degrees,
                                                     control_points,
                                                     weights);
#endif
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateRationalBezier. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
#ifdef SPLINEPY_MORE
  case 4:
    switch (dim) {
    case 1:
      return std::make_shared<RationalBezier<4, 1>>(degrees,
                                                    control_points,
                                                    weights);
    case 2:
      return std::make_shared<RationalBezier<4, 2>>(degrees,
                                                    control_points,
                                                    weights);
    case 3:
      return std::make_shared<RationalBezier<4, 3>>(degrees,
                                                    control_points,
                                                    weights);
    case 4:
      return std::make_shared<RationalBezier<4, 4>>(degrees,
                                                    control_points,
                                                    weights);
    case 5:
      return std::make_shared<RationalBezier<4, 5>>(degrees,
                                                    control_points,
                                                    weights);
    case 6:
      return std::make_shared<RationalBezier<4, 6>>(degrees,
                                                    control_points,
                                                    weights);
    case 7:
      return std::make_shared<RationalBezier<4, 7>>(degrees,
                                                    control_points,
                                                    weights);
    case 8:
      return std::make_shared<RationalBezier<4, 8>>(degrees,
                                                    control_points,
                                                    weights);
    case 9:
      return std::make_shared<RationalBezier<4, 9>>(degrees,
                                                    control_points,
                                                    weights);
    case 10:
      return std::make_shared<RationalBezier<4, 10>>(degrees,
                                                     control_points,
                                                     weights);
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateRationalBezier. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
  case 5:
    switch (dim) {
    case 1:
      return std::make_shared<RationalBezier<5, 1>>(degrees,
                                                    control_points,
                                                    weights);
    case 2:
      return std::make_shared<RationalBezier<5, 2>>(degrees,
                                                    control_points,
                                                    weights);
    case 3:
      return std::make_shared<RationalBezier<5, 3>>(degrees,
                                                    control_points,
                                                    weights);
    case 4:
      return std::make_shared<RationalBezier<5, 4>>(degrees,
                                                    control_points,
                                                    weights);
    case 5:
      return std::make_shared<RationalBezier<5, 5>>(degrees,
                                                    control_points,
                                                    weights);
    case 6:
      return std::make_shared<RationalBezier<5, 6>>(degrees,
                                                    control_points,
                                                    weights);
    case 7:
      return std::make_shared<RationalBezier<5, 7>>(degrees,
                                                    control_points,
                                                    weights);
    case 8:
      return std::make_shared<RationalBezier<5, 8>>(degrees,
                                                    control_points,
                                                    weights);
    case 9:
      return std::make_shared<RationalBezier<5, 9>>(degrees,
                                                    control_points,
                                                    weights);
    case 10:
      return std::make_shared<RationalBezier<5, 10>>(degrees,
                                                     control_points,
                                                     weights);
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateRationalBezier. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
  case 6:
    switch (dim) {
    case 1:
      return std::make_shared<RationalBezier<6, 1>>(degrees,
                                                    control_points,
                                                    weights);
    case 2:
      return std::make_shared<RationalBezier<6, 2>>(degrees,
                                                    control_points,
                                                    weights);
    case 3:
      return std::make_shared<RationalBezier<6, 3>>(degrees,
                                                    control_points,
                                                    weights);
    case 4:
      return std::make_shared<RationalBezier<6, 4>>(degrees,
                                                    control_points,
                                                    weights);
    case 5:
      return std::make_shared<RationalBezier<6, 5>>(degrees,
                                                    control_points,
                                                    weights);
    case 6:
      return std::make_shared<RationalBezier<6, 6>>(degrees,
                                                    control_points,
                                                    weights);
    case 7:
      return std::make_shared<RationalBezier<6, 7>>(degrees,
                                                    control_points,
                                                    weights);
    case 8:
      return std::make_shared<RationalBezier<6, 8>>(degrees,
                                                    control_points,
                                                    weights);
    case 9:
      return std::make_shared<RationalBezier<6, 9>>(degrees,
                                                    control_points,
                                                    weights);
    case 10:
      return std::make_shared<RationalBezier<6, 10>>(degrees,
                                                     control_points,
                                                     weights);
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateRationalBezier. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
  case 7:
    switch (dim) {
    case 1:
      return std::make_shared<RationalBezier<7, 1>>(degrees,
                                                    control_points,
                                                    weights);
    case 2:
      return std::make_shared<RationalBezier<7, 2>>(degrees,
                                                    control_points,
                                                    weights);
    case 3:
      return std::make_shared<RationalBezier<7, 3>>(degrees,
                                                    control_points,
                                                    weights);
    case 4:
      return std::make_shared<RationalBezier<7, 4>>(degrees,
                                                    control_points,
                                                    weights);
    case 5:
      return std::make_shared<RationalBezier<7, 5>>(degrees,
                                                    control_points,
                                                    weights);
    case 6:
      return std::make_shared<RationalBezier<7, 6>>(degrees,
                                                    control_points,
                                                    weights);
    case 7:
      return std::make_shared<RationalBezier<7, 7>>(degrees,
                                                    control_points,
                                                    weights);
    case 8:
      return std::make_shared<RationalBezier<7, 8>>(degrees,
                                                    control_points,
                                                    weights);
    case 9:
      return std::make_shared<RationalBezier<7, 9>>(degrees,
                                                    control_points,
                                                    weights);
    case 10:
      return std::make_shared<RationalBezier<7, 10>>(degrees,
                                                     control_points,
                                                     weights);
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateRationalBezier. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
  case 8:
    switch (dim) {
    case 1:
      return std::make_shared<RationalBezier<8, 1>>(degrees,
                                                    control_points,
                                                    weights);
    case 2:
      return std::make_shared<RationalBezier<8, 2>>(degrees,
                                                    control_points,
                                                    weights);
    case 3:
      return std::make_shared<RationalBezier<8, 3>>(degrees,
                                                    control_points,
                                                    weights);
    case 4:
      return std::make_shared<RationalBezier<8, 4>>(degrees,
                                                    control_points,
                                                    weights);
    case 5:
      return std::make_shared<RationalBezier<8, 5>>(degrees,
                                                    control_points,
                                                    weights);
    case 6:
      return std::make_shared<RationalBezier<8, 6>>(degrees,
                                                    control_points,
                                                    weights);
    case 7:
      return std::make_shared<RationalBezier<8, 7>>(degrees,
                                                    control_points,
                                                    weights);
    case 8:
      return std::make_shared<RationalBezier<8, 8>>(degrees,
                                                    control_points,
                                                    weights);
    case 9:
      return std::make_shared<RationalBezier<8, 9>>(degrees,
                                                    control_points,
                                                    weights);
    case 10:
      return std::make_shared<RationalBezier<8, 10>>(degrees,
                                                     control_points,
                                                     weights);
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateRationalBezier. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
  case 9:
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
          "Something went wrong during CreateRationalBezier. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
  case 10:
    switch (dim) {
    case 1:
      return std::make_shared<RationalBezier<10, 1>>(degrees,
                                                     control_points,
                                                     weights);
    case 2:
      return std::make_shared<RationalBezier<10, 2>>(degrees,
                                                     control_points,
                                                     weights);
    case 3:
      return std::make_shared<RationalBezier<10, 3>>(degrees,
                                                     control_points,
                                                     weights);
    case 4:
      return std::make_shared<RationalBezier<10, 4>>(degrees,
                                                     control_points,
                                                     weights);
    case 5:
      return std::make_shared<RationalBezier<10, 5>>(degrees,
                                                     control_points,
                                                     weights);
    case 6:
      return std::make_shared<RationalBezier<10, 6>>(degrees,
                                                     control_points,
                                                     weights);
    case 7:
      return std::make_shared<RationalBezier<10, 7>>(degrees,
                                                     control_points,
                                                     weights);
    case 8:
      return std::make_shared<RationalBezier<10, 8>>(degrees,
                                                     control_points,
                                                     weights);
    case 9:
      return std::make_shared<RationalBezier<10, 9>>(degrees,
                                                     control_points,
                                                     weights);
    case 10:
      return std::make_shared<RationalBezier<10, 10>>(degrees,
                                                      control_points,
                                                      weights);
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateRationalBezier. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
#endif
  default:
    splinepy::utils::PrintAndThrowError(
        "Something went wrong during CreateRationalBezier. Please help us by "
        "writing "
        "an issue about this case at [ github.com/tataratat/splinepy ]");
  }
  splinepy::utils::PrintAndThrowError(
      "Something went very wrong during CreateRationalBezier. Please help us "
      "by writing "
      "an issue about this case at [ github.com/tataratat/splinepy ]");
  // make compiler happy
  return std::shared_ptr<SplinepyBase>{};
}

} // namespace splinepy::splines
