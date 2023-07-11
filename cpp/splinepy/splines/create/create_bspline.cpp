#include <splinepy/splines/create/create_bspline.hpp>

namespace splinepy::splines::create {

/// dynamic creation of templated Bezier
std::shared_ptr<splinepy::splines::SplinepyBase>
CreateBSpline(const int para_dim,
              const int dim,
              const int* degrees,
              const std::vector<std::vector<double>>* knot_vectors,
              const double* control_points) {
  switch (para_dim) {
  case 1:
    return std::make_shared<BSpline<1>>(degrees,
                                        *knot_vectors,
                                        control_points,
                                        dim);
  case 2:
    return std::make_shared<BSpline<2>>(degrees,
                                        *knot_vectors,
                                        control_points,
                                        dim);
  case 3:
    return std::make_shared<BSpline<3>>(degrees,
                                        *knot_vectors,
                                        control_points,
                                        dim);
#ifdef SPLINEPY_MORE
  case 4:
    return std::make_shared<BSpline<4>>(degrees,
                                        *knot_vectors,
                                        control_points,
                                        dim);
  case 5:
    return std::make_shared<BSpline<5>>(degrees,
                                        *knot_vectors,
                                        control_points,
                                        dim);
  case 6:
    return std::make_shared<BSpline<6>>(degrees,
                                        *knot_vectors,
                                        control_points,
                                        dim);
  case 7:
    return std::make_shared<BSpline<7>>(degrees,
                                        *knot_vectors,
                                        control_points,
                                        dim);
  case 8:
    return std::make_shared<BSpline<8>>(degrees,
                                        *knot_vectors,
                                        control_points,
                                        dim);
  case 9:
    return std::make_shared<BSpline<9>>(degrees,
                                        *knot_vectors,
                                        control_points,
                                        dim);
  case 10:
    return std::make_shared<BSpline<10>>(degrees,
                                         *knot_vectors,
                                         control_points,
                                         dim);
#endif
  default:
    splinepy::utils::PrintAndThrowError(
        "Something went wrong during CreateBSpline. Please help us by writing "
        "an issue about this case at [ github.com/tataratat/splinepy ]");
  }
  splinepy::utils::PrintAndThrowError(
      "Something went very wrong during CreateBSpline. Please help us by "
      "writing "
      "an issue about this case at [ github.com/tataratat/splinepy ]");
  // make compiler happy
  return std::shared_ptr<SplinepyBase>{};
}

} // namespace splinepy::splines::create
