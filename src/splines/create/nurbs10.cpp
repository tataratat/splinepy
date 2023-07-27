#ifdef SPLINEPY_MORE

#include <splinepy/splines/create/create_nurbs.hpp>

namespace splinepy::splines::create {

/// dynamic creation of templated Bezier
std::shared_ptr<splinepy::splines::SplinepyBase>
CreateNurbs10(const int dim,
              const int* degrees,
              const std::vector<std::vector<double>>* knot_vectors,
              const double* control_points,
              const double* weights) {
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
