#include <splinepy/splines/create/create_bezier.hpp>
#include <splinepy/splines/create/create_bspline.hpp>
#include <splinepy/splines/create/create_nurbs.hpp>
#include <splinepy/splines/create/create_rational_bezier.hpp>
#include <splinepy/splines/splinepy_base.hpp>

namespace splinepy::splines {

/// dynamic creation of templated Bezier
std::shared_ptr<SplinepyBase>
splinepy::splines::SplinepyBase::SplinepyCreateBezier(
    const int para_dim,
    const int dim,
    const int* degrees,
    const double* control_points) {
  switch (para_dim) {
  case 1:
    return splinepy::splines::create::CreateBezier1(dim,
                                                    degrees,
                                                    control_points);
  case 2:
    return splinepy::splines::create::CreateBezier2(dim,
                                                    degrees,
                                                    control_points);
  case 3:
    return splinepy::splines::create::CreateBezier3(dim,
                                                    degrees,
                                                    control_points);
#ifdef SPLINEPY_MORE

  case 4:
    return splinepy::splines::create::CreateBezier4(dim,
                                                    degrees,
                                                    control_points);
  case 5:
    return splinepy::splines::create::CreateBezier5(dim,
                                                    degrees,
                                                    control_points);
  case 6:
    return splinepy::splines::create::CreateBezier6(dim,
                                                    degrees,
                                                    control_points);
  case 7:
    return splinepy::splines::create::CreateBezier7(dim,
                                                    degrees,
                                                    control_points);
  case 8:
    return splinepy::splines::create::CreateBezier8(dim,
                                                    degrees,
                                                    control_points);
  case 9:
    return splinepy::splines::create::CreateBezier9(dim,
                                                    degrees,
                                                    control_points);
  case 10:
    return splinepy::splines::create::CreateBezier10(dim,
                                                     degrees,
                                                     control_points);
#endif

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
    return splinepy::splines::create::CreateRationalBezier1(dim,
                                                            degrees,
                                                            control_points,
                                                            weights);
  case 2:
    return splinepy::splines::create::CreateRationalBezier2(dim,
                                                            degrees,
                                                            control_points,
                                                            weights);
  case 3:
    return splinepy::splines::create::CreateRationalBezier3(dim,
                                                            degrees,
                                                            control_points,
                                                            weights);
#ifdef SPLINEPY_MORE

  case 4:
    return splinepy::splines::create::CreateRationalBezier4(dim,
                                                            degrees,
                                                            control_points,
                                                            weights);
  case 5:
    return splinepy::splines::create::CreateRationalBezier5(dim,
                                                            degrees,
                                                            control_points,
                                                            weights);
  case 6:
    return splinepy::splines::create::CreateRationalBezier6(dim,
                                                            degrees,
                                                            control_points,
                                                            weights);
  case 7:
    return splinepy::splines::create::CreateRationalBezier7(dim,
                                                            degrees,
                                                            control_points,
                                                            weights);
  case 8:
    return splinepy::splines::create::CreateRationalBezier8(dim,
                                                            degrees,
                                                            control_points,
                                                            weights);
  case 9:
    return splinepy::splines::create::CreateRationalBezier9(dim,
                                                            degrees,
                                                            control_points,
                                                            weights);
  case 10:
    return splinepy::splines::create::CreateRationalBezier10(dim,
                                                             degrees,
                                                             control_points,
                                                             weights);
#endif

  default:
    splinepy::utils::PrintAndThrowError(
        "Something went wrong during CreateRationalBezier. Please help us by "
        "writing "
        "an issue about this case at [ github.com/tataratat/splinepy ]");
    break;
  }
  splinepy::utils::PrintAndThrowError(
      "Something went very wrong during CreateRationalBezier. Please help us "
      "by "
      "writing "
      "an issue about this case at [ github.com/tataratat/splinepy ]");
  // make compiler happy
  return std::shared_ptr<SplinepyBase>{};
}

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
    return splinepy::splines::create::CreateBSpline1(dim,
                                                     degrees,
                                                     knot_vectors,
                                                     control_points);
  case 2:
    return splinepy::splines::create::CreateBSpline2(dim,
                                                     degrees,
                                                     knot_vectors,
                                                     control_points);
  case 3:
    return splinepy::splines::create::CreateBSpline3(dim,
                                                     degrees,
                                                     knot_vectors,
                                                     control_points);
#ifdef SPLINEPY_MORE

  case 4:
    return splinepy::splines::create::CreateBSpline4(dim,
                                                     degrees,
                                                     knot_vectors,
                                                     control_points);
  case 5:
    return splinepy::splines::create::CreateBSpline5(dim,
                                                     degrees,
                                                     knot_vectors,
                                                     control_points);
  case 6:
    return splinepy::splines::create::CreateBSpline6(dim,
                                                     degrees,
                                                     knot_vectors,
                                                     control_points);
  case 7:
    return splinepy::splines::create::CreateBSpline7(dim,
                                                     degrees,
                                                     knot_vectors,
                                                     control_points);
  case 8:
    return splinepy::splines::create::CreateBSpline8(dim,
                                                     degrees,
                                                     knot_vectors,
                                                     control_points);
  case 9:
    return splinepy::splines::create::CreateBSpline9(dim,
                                                     degrees,
                                                     knot_vectors,
                                                     control_points);
  case 10:
    return splinepy::splines::create::CreateBSpline10(dim,
                                                      degrees,
                                                      knot_vectors,
                                                      control_points);
#endif

  default:
    splinepy::utils::PrintAndThrowError(
        "Something went wrong during CreateBSpline. Please help us by writing "
        "an issue about this case at [ github.com/tataratat/splinepy ]");
    break;
  }
  splinepy::utils::PrintAndThrowError(
      "Something went very wrong during CreateBSpline. Please help us by "
      "writing "
      "an issue about this case at [ github.com/tataratat/splinepy ]");
  // make compiler happy
  return std::shared_ptr<SplinepyBase>{};
}

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
    return splinepy::splines::create::CreateNurbs1(dim,
                                                   degrees,
                                                   knot_vectors,
                                                   control_points,
                                                   weights);
  case 2:
    return splinepy::splines::create::CreateNurbs2(dim,
                                                   degrees,
                                                   knot_vectors,
                                                   control_points,
                                                   weights);
  case 3:
    return splinepy::splines::create::CreateNurbs3(dim,
                                                   degrees,
                                                   knot_vectors,
                                                   control_points,
                                                   weights);
#ifdef SPLINEPY_MORE

  case 4:
    return splinepy::splines::create::CreateNurbs4(dim,
                                                   degrees,
                                                   knot_vectors,
                                                   control_points,
                                                   weights);
  case 5:
    return splinepy::splines::create::CreateNurbs5(dim,
                                                   degrees,
                                                   knot_vectors,
                                                   control_points,
                                                   weights);
  case 6:
    return splinepy::splines::create::CreateNurbs6(dim,
                                                   degrees,
                                                   knot_vectors,
                                                   control_points,
                                                   weights);
  case 7:
    return splinepy::splines::create::CreateNurbs7(dim,
                                                   degrees,
                                                   knot_vectors,
                                                   control_points,
                                                   weights);
  case 8:
    return splinepy::splines::create::CreateNurbs8(dim,
                                                   degrees,
                                                   knot_vectors,
                                                   control_points,
                                                   weights);
  case 9:
    return splinepy::splines::create::CreateNurbs9(dim,
                                                   degrees,
                                                   knot_vectors,
                                                   control_points,
                                                   weights);
  case 10:
    return splinepy::splines::create::CreateNurbs10(dim,
                                                    degrees,
                                                    knot_vectors,
                                                    control_points,
                                                    weights);
#endif

  default:
    splinepy::utils::PrintAndThrowError(
        "Something went wrong during CreateNurbs. Please help us by writing "
        "an issue about this case at [ github.com/tataratat/splinepy ]");
    break;
  }
  splinepy::utils::PrintAndThrowError(
      "Something went very wrong during CreateNurbs. Please help us by "
      "writing "
      "an issue about this case at [ github.com/tataratat/splinepy ]");
  // make compiler happy
  return std::shared_ptr<SplinepyBase>{};
}

} // namespace splinepy::splines
