#include <splinepy/splines/bspline.hpp>
#include <splinepy/splines/splinepy_base.hpp>

namespace splinepy::splines {

/// dynamic creation of templated BSpline
std::shared_ptr<SplinepyBase>
splinepy::splines::SplinepyBase::SplinepyCreateBSpline(
    const int para_dim,
    const int dim,
    const int* degrees,
    const std::vector<std::vector<double>>* knot_vectors,
    const double* control_points) {
  switch (para_dim) {
  case 1:
    switch (dim) {
    case 1:
      return std::make_shared<BSpline<1, 1>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 2:
      return std::make_shared<BSpline<1, 2>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 3:
      return std::make_shared<BSpline<1, 3>>(degrees,
                                             *knot_vectors,
                                             control_points);
#ifdef SPLINEPY_MORE
    case 4:
      return std::make_shared<BSpline<1, 4>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 5:
      return std::make_shared<BSpline<1, 5>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 6:
      return std::make_shared<BSpline<1, 6>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 7:
      return std::make_shared<BSpline<1, 7>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 8:
      return std::make_shared<BSpline<1, 8>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 9:
      return std::make_shared<BSpline<1, 9>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 10:
      return std::make_shared<BSpline<1, 10>>(degrees,
                                              *knot_vectors,
                                              control_points);
#endif
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateBSpline. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
  case 2:
    switch (dim) {
    case 1:
      return std::make_shared<BSpline<2, 1>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 2:
      return std::make_shared<BSpline<2, 2>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 3:
      return std::make_shared<BSpline<2, 3>>(degrees,
                                             *knot_vectors,
                                             control_points);
#ifdef SPLINEPY_MORE
    case 4:
      return std::make_shared<BSpline<2, 4>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 5:
      return std::make_shared<BSpline<2, 5>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 6:
      return std::make_shared<BSpline<2, 6>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 7:
      return std::make_shared<BSpline<2, 7>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 8:
      return std::make_shared<BSpline<2, 8>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 9:
      return std::make_shared<BSpline<2, 9>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 10:
      return std::make_shared<BSpline<2, 10>>(degrees,
                                              *knot_vectors,
                                              control_points);
#endif
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateBSpline. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
  case 3:
    switch (dim) {
    case 1:
      return std::make_shared<BSpline<3, 1>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 2:
      return std::make_shared<BSpline<3, 2>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 3:
      return std::make_shared<BSpline<3, 3>>(degrees,
                                             *knot_vectors,
                                             control_points);
#ifdef SPLINEPY_MORE
    case 4:
      return std::make_shared<BSpline<3, 4>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 5:
      return std::make_shared<BSpline<3, 5>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 6:
      return std::make_shared<BSpline<3, 6>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 7:
      return std::make_shared<BSpline<3, 7>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 8:
      return std::make_shared<BSpline<3, 8>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 9:
      return std::make_shared<BSpline<3, 9>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 10:
      return std::make_shared<BSpline<3, 10>>(degrees,
                                              *knot_vectors,
                                              control_points);
#endif
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateBSpline. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
#ifdef SPLINEPY_MORE
  case 4:
    switch (dim) {
    case 1:
      return std::make_shared<BSpline<4, 1>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 2:
      return std::make_shared<BSpline<4, 2>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 3:
      return std::make_shared<BSpline<4, 3>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 4:
      return std::make_shared<BSpline<4, 4>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 5:
      return std::make_shared<BSpline<4, 5>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 6:
      return std::make_shared<BSpline<4, 6>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 7:
      return std::make_shared<BSpline<4, 7>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 8:
      return std::make_shared<BSpline<4, 8>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 9:
      return std::make_shared<BSpline<4, 9>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 10:
      return std::make_shared<BSpline<4, 10>>(degrees,
                                              *knot_vectors,
                                              control_points);
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateBSpline. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
  case 5:
    switch (dim) {
    case 1:
      return std::make_shared<BSpline<5, 1>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 2:
      return std::make_shared<BSpline<5, 2>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 3:
      return std::make_shared<BSpline<5, 3>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 4:
      return std::make_shared<BSpline<5, 4>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 5:
      return std::make_shared<BSpline<5, 5>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 6:
      return std::make_shared<BSpline<5, 6>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 7:
      return std::make_shared<BSpline<5, 7>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 8:
      return std::make_shared<BSpline<5, 8>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 9:
      return std::make_shared<BSpline<5, 9>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 10:
      return std::make_shared<BSpline<5, 10>>(degrees,
                                              *knot_vectors,
                                              control_points);
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateBSpline. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
  case 6:
    switch (dim) {
    case 1:
      return std::make_shared<BSpline<6, 1>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 2:
      return std::make_shared<BSpline<6, 2>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 3:
      return std::make_shared<BSpline<6, 3>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 4:
      return std::make_shared<BSpline<6, 4>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 5:
      return std::make_shared<BSpline<6, 5>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 6:
      return std::make_shared<BSpline<6, 6>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 7:
      return std::make_shared<BSpline<6, 7>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 8:
      return std::make_shared<BSpline<6, 8>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 9:
      return std::make_shared<BSpline<6, 9>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 10:
      return std::make_shared<BSpline<6, 10>>(degrees,
                                              *knot_vectors,
                                              control_points);
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateBSpline. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
  case 7:
    switch (dim) {
    case 1:
      return std::make_shared<BSpline<7, 1>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 2:
      return std::make_shared<BSpline<7, 2>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 3:
      return std::make_shared<BSpline<7, 3>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 4:
      return std::make_shared<BSpline<7, 4>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 5:
      return std::make_shared<BSpline<7, 5>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 6:
      return std::make_shared<BSpline<7, 6>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 7:
      return std::make_shared<BSpline<7, 7>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 8:
      return std::make_shared<BSpline<7, 8>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 9:
      return std::make_shared<BSpline<7, 9>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 10:
      return std::make_shared<BSpline<7, 10>>(degrees,
                                              *knot_vectors,
                                              control_points);
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateBSpline. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
  case 8:
    switch (dim) {
    case 1:
      return std::make_shared<BSpline<8, 1>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 2:
      return std::make_shared<BSpline<8, 2>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 3:
      return std::make_shared<BSpline<8, 3>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 4:
      return std::make_shared<BSpline<8, 4>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 5:
      return std::make_shared<BSpline<8, 5>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 6:
      return std::make_shared<BSpline<8, 6>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 7:
      return std::make_shared<BSpline<8, 7>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 8:
      return std::make_shared<BSpline<8, 8>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 9:
      return std::make_shared<BSpline<8, 9>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 10:
      return std::make_shared<BSpline<8, 10>>(degrees,
                                              *knot_vectors,
                                              control_points);
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateBSpline. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
  case 9:
    switch (dim) {
    case 1:
      return std::make_shared<BSpline<9, 1>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 2:
      return std::make_shared<BSpline<9, 2>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 3:
      return std::make_shared<BSpline<9, 3>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 4:
      return std::make_shared<BSpline<9, 4>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 5:
      return std::make_shared<BSpline<9, 5>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 6:
      return std::make_shared<BSpline<9, 6>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 7:
      return std::make_shared<BSpline<9, 7>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 8:
      return std::make_shared<BSpline<9, 8>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 9:
      return std::make_shared<BSpline<9, 9>>(degrees,
                                             *knot_vectors,
                                             control_points);
    case 10:
      return std::make_shared<BSpline<9, 10>>(degrees,
                                              *knot_vectors,
                                              control_points);
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateBSpline. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
  case 10:
    switch (dim) {
    case 1:
      return std::make_shared<BSpline<10, 1>>(degrees,
                                              *knot_vectors,
                                              control_points);
    case 2:
      return std::make_shared<BSpline<10, 2>>(degrees,
                                              *knot_vectors,
                                              control_points);
    case 3:
      return std::make_shared<BSpline<10, 3>>(degrees,
                                              *knot_vectors,
                                              control_points);
    case 4:
      return std::make_shared<BSpline<10, 4>>(degrees,
                                              *knot_vectors,
                                              control_points);
    case 5:
      return std::make_shared<BSpline<10, 5>>(degrees,
                                              *knot_vectors,
                                              control_points);
    case 6:
      return std::make_shared<BSpline<10, 6>>(degrees,
                                              *knot_vectors,
                                              control_points);
    case 7:
      return std::make_shared<BSpline<10, 7>>(degrees,
                                              *knot_vectors,
                                              control_points);
    case 8:
      return std::make_shared<BSpline<10, 8>>(degrees,
                                              *knot_vectors,
                                              control_points);
    case 9:
      return std::make_shared<BSpline<10, 9>>(degrees,
                                              *knot_vectors,
                                              control_points);
    case 10:
      return std::make_shared<BSpline<10, 10>>(degrees,
                                               *knot_vectors,
                                               control_points);
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateBSpline. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
#endif
  default:
    splinepy::utils::PrintAndThrowError(
        "Something went wrong during CreateBSpline. Please help us by writing "
        "an issue about this case at [ github.com/tataratat/splinepy ]");
  }
  // make compiler happy
  return std::shared_ptr<SplinepyBase>{};
}

} // namespace splinepy::splines
