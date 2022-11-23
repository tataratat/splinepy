#include <splinepy/splines/nurbs.hpp>
#include <splinepy/splines/splinepy_base.hpp>

namespace splinepy::splines {

/// dynamic creation of templated Nurbs
std::shared_ptr<SplinepyBase>
splinepy::splines::SplinepyBase::SplinepyCreateNurbs(
    const int para_dim,
    const int dim,
    const int* degrees,
    const std::vector<std::vector<double>>* knot_vectors,
    const double* control_points,
    const double* weights) {

  switch (para_dim) {
  case 1:
    switch (dim) {
    case 1:
      return std::make_shared<Nurbs<1, 1>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 2:
      return std::make_shared<Nurbs<1, 2>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 3:
      return std::make_shared<Nurbs<1, 3>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
#ifdef SPLINEPY_MORE
    case 4:
      return std::make_shared<Nurbs<1, 4>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 5:
      return std::make_shared<Nurbs<1, 5>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 6:
      return std::make_shared<Nurbs<1, 6>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 7:
      return std::make_shared<Nurbs<1, 7>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 8:
      return std::make_shared<Nurbs<1, 8>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 9:
      return std::make_shared<Nurbs<1, 9>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 10:
      return std::make_shared<Nurbs<1, 10>>(degrees,
                                            *knot_vectors,
                                            control_points,
                                            weights);
#endif
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateNurbs. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
  case 2:
    switch (dim) {
    case 1:
      return std::make_shared<Nurbs<2, 1>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 2:
      return std::make_shared<Nurbs<2, 2>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 3:
      return std::make_shared<Nurbs<2, 3>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
#ifdef SPLINEPY_MORE
    case 4:
      return std::make_shared<Nurbs<2, 4>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 5:
      return std::make_shared<Nurbs<2, 5>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 6:
      return std::make_shared<Nurbs<2, 6>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 7:
      return std::make_shared<Nurbs<2, 7>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 8:
      return std::make_shared<Nurbs<2, 8>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 9:
      return std::make_shared<Nurbs<2, 9>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 10:
      return std::make_shared<Nurbs<2, 10>>(degrees,
                                            *knot_vectors,
                                            control_points,
                                            weights);
#endif
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateNurbs. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
  case 3:
    switch (dim) {
    case 1:
      return std::make_shared<Nurbs<3, 1>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 2:
      return std::make_shared<Nurbs<3, 2>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 3:
      return std::make_shared<Nurbs<3, 3>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
#ifdef SPLINEPY_MORE
    case 4:
      return std::make_shared<Nurbs<3, 4>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 5:
      return std::make_shared<Nurbs<3, 5>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 6:
      return std::make_shared<Nurbs<3, 6>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 7:
      return std::make_shared<Nurbs<3, 7>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 8:
      return std::make_shared<Nurbs<3, 8>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 9:
      return std::make_shared<Nurbs<3, 9>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 10:
      return std::make_shared<Nurbs<3, 10>>(degrees,
                                            *knot_vectors,
                                            control_points,
                                            weights);
#endif
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateNurbs. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
#ifdef SPLINEPY_MORE
  case 4:
    switch (dim) {
    case 1:
      return std::make_shared<Nurbs<4, 1>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 2:
      return std::make_shared<Nurbs<4, 2>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 3:
      return std::make_shared<Nurbs<4, 3>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 4:
      return std::make_shared<Nurbs<4, 4>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 5:
      return std::make_shared<Nurbs<4, 5>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 6:
      return std::make_shared<Nurbs<4, 6>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 7:
      return std::make_shared<Nurbs<4, 7>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 8:
      return std::make_shared<Nurbs<4, 8>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 9:
      return std::make_shared<Nurbs<4, 9>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 10:
      return std::make_shared<Nurbs<4, 10>>(degrees,
                                            *knot_vectors,
                                            control_points,
                                            weights);
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateNurbs. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
  case 5:
    switch (dim) {
    case 1:
      return std::make_shared<Nurbs<5, 1>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 2:
      return std::make_shared<Nurbs<5, 2>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 3:
      return std::make_shared<Nurbs<5, 3>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 4:
      return std::make_shared<Nurbs<5, 4>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 5:
      return std::make_shared<Nurbs<5, 5>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 6:
      return std::make_shared<Nurbs<5, 6>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 7:
      return std::make_shared<Nurbs<5, 7>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 8:
      return std::make_shared<Nurbs<5, 8>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 9:
      return std::make_shared<Nurbs<5, 9>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 10:
      return std::make_shared<Nurbs<5, 10>>(degrees,
                                            *knot_vectors,
                                            control_points,
                                            weights);
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateNurbs. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
  case 6:
    switch (dim) {
    case 1:
      return std::make_shared<Nurbs<6, 1>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 2:
      return std::make_shared<Nurbs<6, 2>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 3:
      return std::make_shared<Nurbs<6, 3>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 4:
      return std::make_shared<Nurbs<6, 4>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 5:
      return std::make_shared<Nurbs<6, 5>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 6:
      return std::make_shared<Nurbs<6, 6>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 7:
      return std::make_shared<Nurbs<6, 7>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 8:
      return std::make_shared<Nurbs<6, 8>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 9:
      return std::make_shared<Nurbs<6, 9>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 10:
      return std::make_shared<Nurbs<6, 10>>(degrees,
                                            *knot_vectors,
                                            control_points,
                                            weights);
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateNurbs. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
  case 7:
    switch (dim) {
    case 1:
      return std::make_shared<Nurbs<7, 1>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 2:
      return std::make_shared<Nurbs<7, 2>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 3:
      return std::make_shared<Nurbs<7, 3>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 4:
      return std::make_shared<Nurbs<7, 4>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 5:
      return std::make_shared<Nurbs<7, 5>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 6:
      return std::make_shared<Nurbs<7, 6>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 7:
      return std::make_shared<Nurbs<7, 7>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 8:
      return std::make_shared<Nurbs<7, 8>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 9:
      return std::make_shared<Nurbs<7, 9>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 10:
      return std::make_shared<Nurbs<7, 10>>(degrees,
                                            *knot_vectors,
                                            control_points,
                                            weights);
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateNurbs. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
  case 8:
    switch (dim) {
    case 1:
      return std::make_shared<Nurbs<8, 1>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 2:
      return std::make_shared<Nurbs<8, 2>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 3:
      return std::make_shared<Nurbs<8, 3>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 4:
      return std::make_shared<Nurbs<8, 4>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 5:
      return std::make_shared<Nurbs<8, 5>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 6:
      return std::make_shared<Nurbs<8, 6>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 7:
      return std::make_shared<Nurbs<8, 7>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 8:
      return std::make_shared<Nurbs<8, 8>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 9:
      return std::make_shared<Nurbs<8, 9>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 10:
      return std::make_shared<Nurbs<8, 10>>(degrees,
                                            *knot_vectors,
                                            control_points,
                                            weights);
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateNurbs. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
  case 9:
    switch (dim) {
    case 1:
      return std::make_shared<Nurbs<9, 1>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 2:
      return std::make_shared<Nurbs<9, 2>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 3:
      return std::make_shared<Nurbs<9, 3>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 4:
      return std::make_shared<Nurbs<9, 4>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 5:
      return std::make_shared<Nurbs<9, 5>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 6:
      return std::make_shared<Nurbs<9, 6>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 7:
      return std::make_shared<Nurbs<9, 7>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 8:
      return std::make_shared<Nurbs<9, 8>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 9:
      return std::make_shared<Nurbs<9, 9>>(degrees,
                                           *knot_vectors,
                                           control_points,
                                           weights);
    case 10:
      return std::make_shared<Nurbs<9, 10>>(degrees,
                                            *knot_vectors,
                                            control_points,
                                            weights);
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateNurbs. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
  case 10:
    switch (dim) {
    case 1:
      return std::make_shared<Nurbs<10, 1>>(degrees,
                                            *knot_vectors,
                                            control_points,
                                            weights);
    case 2:
      return std::make_shared<Nurbs<10, 2>>(degrees,
                                            *knot_vectors,
                                            control_points,
                                            weights);
    case 3:
      return std::make_shared<Nurbs<10, 3>>(degrees,
                                            *knot_vectors,
                                            control_points,
                                            weights);
    case 4:
      return std::make_shared<Nurbs<10, 4>>(degrees,
                                            *knot_vectors,
                                            control_points,
                                            weights);
    case 5:
      return std::make_shared<Nurbs<10, 5>>(degrees,
                                            *knot_vectors,
                                            control_points,
                                            weights);
    case 6:
      return std::make_shared<Nurbs<10, 6>>(degrees,
                                            *knot_vectors,
                                            control_points,
                                            weights);
    case 7:
      return std::make_shared<Nurbs<10, 7>>(degrees,
                                            *knot_vectors,
                                            control_points,
                                            weights);
    case 8:
      return std::make_shared<Nurbs<10, 8>>(degrees,
                                            *knot_vectors,
                                            control_points,
                                            weights);
    case 9:
      return std::make_shared<Nurbs<10, 9>>(degrees,
                                            *knot_vectors,
                                            control_points,
                                            weights);
    case 10:
      return std::make_shared<Nurbs<10, 10>>(degrees,
                                             *knot_vectors,
                                             control_points,
                                             weights);
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during CreateNurbs. Please help us by "
          "writing "
          "an issue about this case at [ github.com/tataratat/splinepy ]");
    }
    break;
#endif
  default:
    splinepy::utils::PrintAndThrowError(
        "Something went wrong during CreateNurbs. Please help us by "
        "writing "
        "an issue about this case at [ github.com/tataratat/splinepy ]");
  }
  splinepy::utils::PrintAndThrowError(
      "Something went very wrong during CreateNurbs. Please help us "
      "by writing "
      "an issue about this case at [ github.com/tataratat/splinepy ]");
  // make compiler happy
  return std::shared_ptr<SplinepyBase>{};
}

} // namespace splinepy::splines
