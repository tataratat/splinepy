#include <memory>
#include <vector>

#include <splinepy/splines/bspline.hpp>
#include <splinepy/splines/create/create_bezier.hpp>
#include <splinepy/splines/create/create_rational_bezier.hpp>
#include <splinepy/splines/nurbs.hpp>
#include <splinepy/splines/splinepy_base.hpp>
#include <splinepy/utils/print.hpp>

namespace splinepy::splines {

std::shared_ptr<SplinepyBase> SplinepyBase::SplinepyCreate(
    const int para_dim,
    const int dim,
    const int* degrees,
    const std::vector<std::vector<double>>* knot_vectors,
    double* control_points,
    double* weights) {
  if (!degrees || !control_points) {
    splinepy::utils::PrintAndThrowError(
        "Not Enough information to create any spline.");
  }

  if (!knot_vectors) {
    if (!weights) {
      return SplinepyCreateBezier(para_dim, dim, degrees, control_points);
    } else {
      return SplinepyCreateRationalBezier(para_dim,
                                          dim,
                                          degrees,
                                          control_points,
                                          weights);
    }
  } else {
    if (!weights) {
      return SplinepyCreateBSpline(para_dim,
                                   dim,
                                   degrees,
                                   knot_vectors,
                                   control_points);
    } else {
      return SplinepyCreateNurbs(para_dim,
                                 dim,
                                 degrees,
                                 knot_vectors,
                                 control_points,
                                 weights);
    }
  }
}

std::shared_ptr<SplinepyBase>
SplinepyBase::SplinepyCreateBezier(const int para_dim,
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

std::shared_ptr<SplinepyBase>
SplinepyBase::SplinepyCreateRationalBezier(const int para_dim,
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

std::shared_ptr<SplinepyBase> SplinepyBase::SplinepyCreateBSpline(
    const int para_dim,
    const int dim,
    const int* degrees,
    const std::vector<std::vector<double>>* knot_vectors,
    double* control_points) {

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

std::shared_ptr<SplinepyBase> SplinepyBase::SplinepyCreateNurbs(
    const int para_dim,
    const int dim,
    const int* degrees,
    const std::vector<std::vector<double>>* knot_vectors,
    double* control_points,
    double* weights) {
  switch (para_dim) {
  case 1:
    return std::make_shared<Nurbs<1>>(degrees,
                                      *knot_vectors,
                                      control_points,
                                      weights,
                                      dim);
  case 2:
    return std::make_shared<Nurbs<2>>(degrees,
                                      *knot_vectors,
                                      control_points,
                                      weights,
                                      dim);
  case 3:
    return std::make_shared<Nurbs<3>>(degrees,
                                      *knot_vectors,
                                      control_points,
                                      weights,
                                      dim);
#ifdef SPLINEPY_MORE
  case 4:
    return std::make_shared<Nurbs<4>>(degrees,
                                      *knot_vectors,
                                      control_points,
                                      weights,
                                      dim);
  case 5:
    return std::make_shared<Nurbs<5>>(degrees,
                                      *knot_vectors,
                                      control_points,
                                      weights,
                                      dim);
  case 6:
    return std::make_shared<Nurbs<6>>(degrees,
                                      *knot_vectors,
                                      control_points,
                                      weights,
                                      dim);
  case 7:
    return std::make_shared<Nurbs<7>>(degrees,
                                      *knot_vectors,
                                      control_points,
                                      weights,
                                      dim);
  case 8:
    return std::make_shared<Nurbs<8>>(degrees,
                                      *knot_vectors,
                                      control_points,
                                      weights,
                                      dim);
  case 9:
    return std::make_shared<Nurbs<9>>(degrees,
                                      *knot_vectors,
                                      control_points,
                                      weights,
                                      dim);
  case 10:
    return std::make_shared<Nurbs<10>>(degrees,
                                       *knot_vectors,
                                       control_points,
                                       weights,
                                       dim);
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

bool SplinepyBase::SplinepySplineNameMatches(const SplinepyBase& a,
                                             const SplinepyBase& b,
                                             const std::string description,
                                             const bool raise) {
  if (a.SplinepySplineName() != b.SplinepySplineName()) {
    if (raise) {
      splinepy::utils::PrintAndThrowError(description,
                                          "Spline name mismatch -"
                                          "Spline0:",
                                          a.SplinepySplineName(),
                                          "/",
                                          "Spline1:",
                                          b.SplinepySplineName());
    }
    return false;
  }
  return true;
}

bool SplinepyBase::SplinepyParaDimMatches(const SplinepyBase& a,
                                          const SplinepyBase& b,
                                          const std::string description,
                                          const bool raise) {
  if (a.SplinepyParaDim() != b.SplinepyParaDim()) {
    if (raise) {
      splinepy::utils::PrintAndThrowError(
          description,
          "Spline parametric dimension mismatch - "
          "Spline0:",
          a.SplinepyParaDim(),
          "/",
          "Spline1:",
          b.SplinepyParaDim());
    }
    return false;
  }
  return true;
}

bool SplinepyBase::SplinepyDimMatches(const SplinepyBase& a,
                                      const SplinepyBase& b,
                                      const std::string description,
                                      const bool raise) {
  if (a.SplinepyDim() != b.SplinepyDim()) {
    if (raise) {
      splinepy::utils::PrintAndThrowError(
          description,
          "Spline parametric dimension mismatch - "
          "Spline0:",
          a.SplinepyDim(),
          "/",
          "Spline1:",
          b.SplinepyDim());
    }
    return false;
  }
  return true;
}

std::shared_ptr<bsplinelib::parameter_spaces::KnotVector>
SplinepyBase::SplinepyKnotVector(const int p_dim) {
  splinepy::utils::PrintAndThrowError("SplinepyKnotVector not implemented for",
                                      SplinepyWhatAmI());
  return nullptr;
}

std::shared_ptr<typename SplinepyBase::ControlPointPointers_>
SplinepyBase::SplinepyControlPointPointers() {
  splinepy::utils::PrintAndThrowError(
      "SplinepyControlPointPointers not implemented for",
      SplinepyWhatAmI());
  return nullptr;
}

std::shared_ptr<typename SplinepyBase::WeightedControlPointPointers_>
SplinepyBase::SplinepyWeightedControlPointPointers() {
  splinepy::utils::PrintAndThrowError(
      "SplinepyWeightedControlPointPointers not implemented for",
      SplinepyWhatAmI());
  return nullptr;
}

std::shared_ptr<typename SplinepyBase::WeightPointers_>
SplinepyBase::SplinepyWeightPointers() {
  splinepy::utils::PrintAndThrowError(
      "SplinepyWeightPointers not implemented for",
      SplinepyWhatAmI());
  return nullptr;
}

void SplinepyBase::SplinepyParametricBounds(double* para_bounds) const {
  splinepy::utils::PrintAndThrowError(
      "SplinepyParametricBounds not implemented for",
      SplinepyWhatAmI());
}

void SplinepyBase::SplinepyControlMeshResolutions(int* control_mesh_res) const {
  splinepy::utils::PrintAndThrowError(
      "SplinepyControlMeshResolutions not implemented for",
      SplinepyWhatAmI());
}

void SplinepyBase::SplinepyGrevilleAbscissae(
    double* greville_abscissae,
    const int& i_para_dim,
    const double& duplicate_tolerance) const {
  splinepy::utils::PrintAndThrowError(
      "SplinepyGrevilleAbscissae not implemented for",
      SplinepyWhatAmI());
}

void SplinepyBase::SplinepyEvaluate(const double* para_coord,
                                    double* evaluated) const {
  splinepy::utils::PrintAndThrowError("SplinepyEvaluate not implemented for",
                                      SplinepyWhatAmI());
}

void SplinepyBase::SplinepyDerivative(const double* para_coord,
                                      const int* orders,
                                      double* derived) const {
  splinepy::utils::PrintAndThrowError("SplinepyDerivative not implemented for",
                                      SplinepyWhatAmI());
}

void SplinepyBase::SplinepyJacobian(const double* para_coord,
                                    double* jacobians) const {
  splinepy::utils::PrintAndThrowError("SplinepyJacobian not implemented for",
                                      SplinepyWhatAmI());
}

void SplinepyBase::SplinepyBasis(const double* para_coord,
                                 double* basis) const {
  splinepy::utils::PrintAndThrowError("SplinepyBasis not implemented for",
                                      SplinepyWhatAmI());
}

void SplinepyBase::SplinepyBasisDerivative(const double* para_coord,
                                           const int* order,
                                           double* basis) const {
  splinepy::utils::PrintAndThrowError(
      "SplinepyBasisDerivative not implemented for",
      SplinepyWhatAmI());
}

void SplinepyBase::SplinepySupport(const double* para_coord,
                                   int* support) const {
  splinepy::utils::PrintAndThrowError("SplinepySupport not implemented for",
                                      SplinepyWhatAmI());
}

void SplinepyBase::SplinepyBasisAndSupport(const double* para_coord,
                                           double* basis,
                                           int* support) const {
  splinepy::utils::PrintAndThrowError(
      "SplinepyBasisAndSupport not implemented for",
      SplinepyWhatAmI());
}

void SplinepyBase::SplinepyBasisDerivativeAndSupport(const double* para_coord,
                                                     const int* orders,
                                                     double* basis,
                                                     int* support) const {
  splinepy::utils::PrintAndThrowError(
      "SplinepyBasisDerivativeAndSupport not implemented for",
      SplinepyWhatAmI());
}

void SplinepyBase::SplinepyPlantNewKdTreeForProximity(const int* resolutions,
                                                      const int& nthreads) {
  splinepy::utils::PrintAndThrowError(
      "SplinepyPlantNewKdTreeForProximity not implemented for",
      SplinepyWhatAmI());
}

void SplinepyBase::SplinepyVerboseProximity(const double* query,
                                            const double& tolerance,
                                            const int& max_iterations,
                                            const bool aggressive_bounds,
                                            double* para_coord,
                                            double* phys_coord,
                                            double* phys_diff,
                                            double& distance,
                                            double& convergence_norm,
                                            double* first_derivatives,
                                            double* second_derivatives) const {
  splinepy::utils::PrintAndThrowError(
      "SplinepyVerboseProximity not implemented for",
      SplinepyWhatAmI());
}

void SplinepyBase::SplinepyElevateDegree(const int& para_dims) {
  splinepy::utils::PrintAndThrowError(
      "SplinepyElevateDegree not implemented for",
      SplinepyWhatAmI());
}

bool SplinepyBase::SplinepyReduceDegree(const int& para_dims,
                                        const double& tolerance) {
  splinepy::utils::PrintAndThrowError(
      "SplinepyReduceDegree not implemented for",
      SplinepyWhatAmI());
  return false;
}

bool SplinepyBase::SplinepyInsertKnot(const int& para_dim, const double& knot) {
  splinepy::utils::PrintAndThrowError("SplinepyInsertKnot not implemented for",
                                      SplinepyWhatAmI());
  return false;
}

bool SplinepyBase::SplinepyRemoveKnot(const int& para_dim,
                                      const double& knot,
                                      const double& tolerance) {
  splinepy::utils::PrintAndThrowError("SplinepyRemoveKnot not implemented for",
                                      SplinepyWhatAmI());
  return false;
}

std::vector<std::vector<int>> SplinepyBase::SplinepyKnotMultiplicities() const {
  splinepy::utils::PrintAndThrowError(
      "SplinepyKnotMultiplicities not implemented for",
      SplinepyWhatAmI());
  return std::vector<std::vector<int>>{};
};

std::shared_ptr<SplinepyBase>
SplinepyBase::SplinepyMultiply(const std::shared_ptr<SplinepyBase>& a) const {
  splinepy::utils::PrintAndThrowError("SplinepyMultiply not implemented for",
                                      SplinepyWhatAmI());
  return std::shared_ptr<SplinepyBase>{};
}

std::shared_ptr<SplinepyBase>
SplinepyBase::SplinepyAdd(const std::shared_ptr<SplinepyBase>& a) const {
  splinepy::utils::PrintAndThrowError("SplinepyAdd not implemented for",
                                      SplinepyWhatAmI());
  return std::shared_ptr<SplinepyBase>{};
}

std::shared_ptr<SplinepyBase> SplinepyBase::SplinepyCompose(
    const std::shared_ptr<SplinepyBase>& inner_function) const {
  splinepy::utils::PrintAndThrowError("SplinepyCompose not implemented for",
                                      SplinepyWhatAmI());
  return std::shared_ptr<SplinepyBase>{};
}

std::vector<std::shared_ptr<SplinepyBase>>
SplinepyBase::SplinepyComposeSensitivities(
    const std::shared_ptr<SplinepyBase>& inner_function) const {
  splinepy::utils::PrintAndThrowError(
      "SplinepyComposeSensitivities not implemented for",
      SplinepyWhatAmI());
  return std::vector<std::shared_ptr<SplinepyBase>>{};
}

std::vector<std::shared_ptr<SplinepyBase>>
SplinepyBase::SplinepySplit(const int& para_dim, const double& location) const {
  splinepy::utils::PrintAndThrowError("SplinepySplit not implemented for",
                                      SplinepyWhatAmI());
  return {std::shared_ptr<SplinepyBase>{}};
}

std::shared_ptr<SplinepyBase>
SplinepyBase::SplinepyDerivativeSpline(const int* orders) const {
  splinepy::utils::PrintAndThrowError(
      "SplinepyDerivativeSpline is not implemented for",
      SplinepyWhatAmI());
  return std::shared_ptr<SplinepyBase>{};
}

std::vector<std::shared_ptr<SplinepyBase>>
SplinepyBase::SplinepyExtractBezierPatches() const {
  splinepy::utils::PrintAndThrowError(
      "SplinepyBezierPatchExtraction is not implemented for",
      SplinepyWhatAmI());
  return {std::shared_ptr<SplinepyBase>{}};
}

std::shared_ptr<SplinepyBase>
SplinepyBase::SplinepyExtractBoundary(const int& boundary_id) {
  splinepy::utils::PrintAndThrowError(
      "SplinepyExtractBoundary is not implemented for",
      SplinepyWhatAmI());
  return {std::shared_ptr<SplinepyBase>{}};
}

std::shared_ptr<SplinepyBase>
SplinepyBase::SplinepyExtractDim(const int& phys_dim) const {
  splinepy::utils::PrintAndThrowError(
      "SplinepyExtractDim is not implemented for",
      SplinepyWhatAmI());
  return {std::shared_ptr<SplinepyBase>{}};
}

std::shared_ptr<SplinepyBase> SplinepyBase::SplinepyCompositionDerivative(
    const std::shared_ptr<SplinepyBase>& inner,
    const std::shared_ptr<SplinepyBase>& inner_derivative) const {
  splinepy::utils::PrintAndThrowError(
      "SplinepyCompositionDerivative is not implemented for",
      SplinepyWhatAmI());
  return {std::shared_ptr<SplinepyBase>{}};
}

std::shared_ptr<SplinepyBase> SplinepyBase::SplinepyDeepCopy() const {
  splinepy::utils::PrintAndThrowError("SplinepyDeepCopy is not implemented for",
                                      SplinepyWhatAmI());
  return {std::shared_ptr<SplinepyBase>{}};
};

} // namespace splinepy::splines
