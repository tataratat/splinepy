#include <array>
#include <type_traits>

#include "splinepy/splines/rational_bezier.hpp"

#include "splinepy/splines/helpers/basis_functions.hpp"
#include "splinepy/splines/helpers/extract.hpp"
#include "splinepy/splines/helpers/properties.hpp"
#include "splinepy/splines/helpers/scalar_type_wrapper.hpp"
#include "splinepy/utils/print.hpp"

namespace splinepy::splines {

template<std::size_t para_dim, std::size_t dim>
typename Bezier<para_dim, dim>::Base_
Bezier<para_dim, dim>::CreateBase(const int* degrees,
                                  const double* control_points) {

  std::array<std::size_t, para_dim> bm_degrees{};
  std::size_t ncps{1};

  // formulate degrees
  for (std::size_t i{}; i < para_dim; ++i) {
    bm_degrees[i] = degrees[i];
    ncps *= degrees[i] + 1;
  }

  // formulate control_points
  std::vector<Coordinate_> bm_control_points(ncps);
  for (std::size_t i{}; i < ncps; ++i) {
    if constexpr (dim > 1) {
      for (std::size_t j = 0; j < dim; j++) {
        bm_control_points[i][j] = control_points[i * dim + j];
      }
    } else {
      bm_control_points[i] = control_points[i];
    }
  }
  return Base_(bm_degrees, bm_control_points);
}
template<std::size_t para_dim, std::size_t dim>
constexpr auto
Bezier<para_dim, dim>::operator()(const ParametricCoordinate_& query) const {
  return Base_::Evaluate(query);
}

template<std::size_t para_dim, std::size_t dim>
constexpr auto
Bezier<para_dim, dim>::operator()(const ParametricCoordinate_& query,
                                  const Derivative_& order) const {
  return Base_::EvaluateDerivative(query, order);
}

template<std::size_t para_dim, std::size_t dim>
constexpr auto Bezier<para_dim, dim>::ElevateDegree(const Dimension_ p_dim) {
  return Base_::OrderElevateAlongParametricDimension(p_dim);
}

template<std::size_t para_dim, std::size_t dim>
int Bezier<para_dim, dim>::SplinepyNumberOfControlPoints() const {
  return static_cast<int>(Base_::control_points.size());
}

template<std::size_t para_dim, std::size_t dim>
int Bezier<para_dim, dim>::SplinepyNumberOfSupports() const {
  return splinepy::splines::helpers::GetNumberOfSupports(*this);
}

template<std::size_t para_dim, std::size_t dim>
void Bezier<para_dim, dim>::SplinepyCurrentProperties(
    int* degrees,
    std::vector<std::vector<double>>* knot_vectors /* untouched */,
    double* control_points,
    double* weights /* untouched */) const {

  // degrees
  if (degrees) {
    for (std::size_t i{}; i < kParaDim; ++i) {
      degrees[i] = static_cast<int>(Base_::GetDegrees()[i]);
    }
  }

  // control_points
  if (control_points) {
    const std::size_t ncps = Base_::control_points.size();
    for (std::size_t i{}; i < ncps; ++i) {
      if constexpr (dim > 1) {
        for (std::size_t j{}; j < kDim; ++j) {
          control_points[i * kDim + j] = Base_::control_points[i][j];
        }
      } else {
        control_points[i] = Base_::control_points[i];
      }
    }
  }
}

template<std::size_t para_dim, std::size_t dim>
std::shared_ptr<typename Bezier<para_dim, dim>::ControlPointPointers_>
Bezier<para_dim, dim>::SplinepyControlPointPointers() {
  if (SplinepyBase_::control_point_pointers_
      && SplinepyBase_::control_point_pointers_->Len()
             == SplinepyNumberOfControlPoints()) {
    return SplinepyBase_::control_point_pointers_;
  }
  auto cpp = std::make_shared<ControlPointPointers_>();
  cpp->dim_ = kDim;
  cpp->coordinate_begins_.reserve(SplinepyNumberOfControlPoints());
  if constexpr (kDim == 1) {
    for (auto& control_point : Base_::control_points) {
      cpp->coordinate_begins_.push_back(&control_point);
    }
  } else {
    for (auto& control_point : Base_::control_points) {
      cpp->coordinate_begins_.push_back(control_point.data());
    }
  }

  SplinepyBase_::control_point_pointers_ = cpp;

  return cpp;
}

template<std::size_t para_dim, std::size_t dim>
void Bezier<para_dim, dim>::SplinepyParametricBounds(
    double* para_bounds) const {
  for (std::size_t i{}; i < kParaDim; ++i) {
    // lower bounds
    para_bounds[i] = 0.;
    // upper bounds
    para_bounds[kParaDim + i] = 1.;
  }
}

template<std::size_t para_dim, std::size_t dim>
void Bezier<para_dim, dim>::SplinepyControlMeshResolutions(
    int* control_mesh_res) const {
  const auto cm_res =
      splinepy::splines::helpers::GetControlMeshResolutions(*this);
  std::copy_n(cm_res.begin(), para_dim, control_mesh_res);
}

template<std::size_t para_dim, std::size_t dim>
void Bezier<para_dim, dim>::SplinepyGrevilleAbscissae(
    double* greville_abscissae,
    const int& i_para_dim,
    const double& duplicate_tolerance) const {
  splinepy::splines::helpers::GetGrevilleAbscissae(*this,
                                                   greville_abscissae,
                                                   i_para_dim,
                                                   duplicate_tolerance);
}
template<std::size_t para_dim, std::size_t dim>
void Bezier<para_dim, dim>::SplinepyEvaluate(const double* para_coord,
                                             double* evaluated) const {
  splinepy::splines::helpers::ScalarTypeEvaluate(*this, para_coord, evaluated);
}

template<std::size_t para_dim, std::size_t dim>
void Bezier<para_dim, dim>::SplinepyDerivative(const double* para_coord,
                                               const int* orders,
                                               double* derived) const {
  splinepy::splines::helpers::ScalarTypeDerivative(*this,
                                                   para_coord,
                                                   orders,
                                                   derived);
}

template<std::size_t para_dim, std::size_t dim>
void Bezier<para_dim, dim>::SplinepyJacobian(const double* para_coord,
                                             double* jacobians) const {
  splinepy::splines::helpers::ScalarTypeJacobian(*this, para_coord, jacobians);
}

template<std::size_t para_dim, std::size_t dim>
void Bezier<para_dim, dim>::SplinepyPlantNewKdTreeForProximity(
    const int* resolutions,
    const int& nthreads) {
  GetProximity().PlantNewKdTree(resolutions, nthreads);
}

template<std::size_t para_dim, std::size_t dim>
void Bezier<para_dim, dim>::SplinepyVerboseProximity(
    const double* query,
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
  GetProximity().VerboseQuery(query,
                              tolerance,
                              max_iterations,
                              aggressive_bounds,
                              para_coord,
                              phys_coord,
                              phys_diff,
                              distance,
                              convergence_norm,
                              first_derivatives,
                              second_derivatives);
}

template<std::size_t para_dim, std::size_t dim>
void Bezier<para_dim, dim>::SplinepyElevateDegree(const int& p_dim) {
  splinepy::splines::helpers::ScalarTypeElevateDegree(*this, p_dim);
}

template<std::size_t para_dim, std::size_t dim>
void Bezier<para_dim, dim>::SplinepyBasis(const double* para_coord,
                                          double* basis) const {
  splinepy::splines::helpers::BezierBasis(*this, para_coord, basis);
}

template<std::size_t para_dim, std::size_t dim>
void Bezier<para_dim, dim>::SplinepyBasisDerivative(const double* para_coord,
                                                    const int* order,
                                                    double* basis_der) const {
  splinepy::splines::helpers::BezierBasisDerivative(*this,
                                                    para_coord,
                                                    order,
                                                    basis_der);
}

template<std::size_t para_dim, std::size_t dim>
void Bezier<para_dim, dim>::SplinepySupport(const double* para_coord,
                                            int* support) const {
  splinepy::splines::helpers::BezierSupport(*this, para_coord, support);
}

template<std::size_t para_dim, std::size_t dim>
void Bezier<para_dim, dim>::SplinepyBasisAndSupport(const double* para_coord,
                                                    double* basis,
                                                    int* support) const {

  SplinepyBasis(para_coord, basis);
  SplinepySupport(para_coord, support);
}

template<std::size_t para_dim, std::size_t dim>

void Bezier<para_dim, dim>::SplinepyBasisDerivativeAndSupport(
    const double* para_coord,
    const int* orders,
    double* basis_der,
    int* support) const {
  SplinepyBasisDerivative(para_coord, orders, basis_der);
  SplinepySupport(para_coord, support);
}

template<std::size_t para_dim, std::size_t dim>
std::shared_ptr<SplinepyBase> Bezier<para_dim, dim>::SplinepyMultiply(
    const std::shared_ptr<SplinepyBase>& a) const {

  SplinepySplineNameMatches(
      *this,
      *a,
      "Spline multiplication requires splines of same type.",
      true);
  SplinepyParaDimMatches(
      *this,
      *a,
      "Spline multiplication requires splines of same parametric dimension.",
      true);
  // check dimension if a is not a scalar spline
  if (a->SplinepyDim() != 1) {
    SplinepyDimMatches(
        *this,
        *a,
        (std::string) "Spline multiplication requires splines of either 1 or "
            + "same physical dimension.",
        true);
  }

  // good to multiply
  if (a->SplinepyDim() == 1) {
    auto true_a = static_cast<Bezier<para_dim, 1>&>(*a);
    return std::make_shared<Bezier<para_dim, dim>>(
        this->Base_::operator*(true_a));
  } else {
    auto true_a = static_cast<Bezier<para_dim, dim>&>(*a);
    return std::make_shared<Bezier<para_dim, 1>>(
        this->Base_::operator*(true_a));
  }
}

template<std::size_t para_dim, std::size_t dim>
std::shared_ptr<SplinepyBase> Bezier<para_dim, dim>::SplinepyAdd(
    const std::shared_ptr<SplinepyBase>& a) const {
  SplinepySplineNameMatches(
      *this,
      *a,
      "Spline addition requires splines of the same type.",
      true);
  SplinepyParaDimMatches(
      *this,
      *a,
      "Spline addition requires splines of the same parametric dimension.",
      true);
  SplinepyDimMatches(
      *this,
      *a,
      "Spline addition requires splines of the same physical dimension",
      true);

  return std::make_shared<Bezier<para_dim, dim>>(
      this->Base_::operator+(static_cast<Bezier<para_dim, dim>&>(*a)));
}

template<std::size_t para_dim, std::size_t dim>
template<bool is_rational, size_t inner_para_dim>
std::vector<std::shared_ptr<SplinepyBase>>
Bezier<para_dim, dim>::ConvertComposeToBase(
    const std::shared_ptr<SplinepyBase>& inner_function) const {
  if constexpr (is_rational) {
    auto super_return = this->ComposeSensitivity(
        static_cast<RationalBezier<inner_para_dim, para_dim>&>(
            *inner_function));

    std::vector<std::shared_ptr<SplinepyBase>> to_return(super_return.size());

    for (size_t i{}; i < super_return.size(); ++i) {
      to_return[i] =
          std::make_shared<RationalBezier<inner_para_dim, 1>>(super_return[i]);
    }
    return to_return;
  } else {
    auto super_return = this->ComposeSensitivity(
        static_cast<Bezier<inner_para_dim, para_dim>&>(*inner_function));

    std::vector<std::shared_ptr<SplinepyBase>> to_return(super_return.size());

    for (size_t i{}; i < super_return.size(); ++i) {
      to_return[i] =
          std::make_shared<Bezier<inner_para_dim, 1>>(super_return[i]);
    }
    return to_return;
  }
}

template<std::size_t para_dim, std::size_t dim>
std::vector<std::shared_ptr<SplinepyBase>>
Bezier<para_dim, dim>::SplinepyComposeSensitivities(
    const std::shared_ptr<SplinepyBase>& inner_function) const {
  // type check
  if (inner_function->SplinepySplineName().find("Bezier")
      == std::string::npos) {
    splinepy::utils::PrintAndThrowError(
        "Bezier composition requires inner function to be a bezier type.",
        "Given inner function -",
        inner_function->SplinepyWhatAmI());
  }

  // composable?
  if (inner_function->SplinepyDim() != this->SplinepyParaDim()) {
    splinepy::utils::PrintAndThrowError(
        "Spline composition requires inner function to have same physical",
        "dimension as outer function's parametric dimension.",
        "Outer Function:",
        this->SplinepyWhatAmI(),
        "/",
        "Inner Function:",
        inner_function->SplinepyWhatAmI());
  }

  // compose - return correct type.
  if (inner_function->SplinepyIsRational()) {
    switch (inner_function->SplinepyParaDim()) {
    case 1:
      return ConvertComposeToBase<true, 1>(inner_function);
    case 2:
      return ConvertComposeToBase<true, 2>(inner_function);
    case 3:
      return ConvertComposeToBase<true, 3>(inner_function);
#ifdef SPLINEPY_MORE
    case 4:
      return ConvertComposeToBase<true, 4>(inner_function);
    case 5:
      return ConvertComposeToBase<true, 5>(inner_function);
    case 6:
      return ConvertComposeToBase<true, 6>(inner_function);
    case 7:
      return ConvertComposeToBase<true, 7>(inner_function);
    case 8:
      return ConvertComposeToBase<true, 8>(inner_function);
    case 9:
      return ConvertComposeToBase<true, 9>(inner_function);
    case 10:
      return ConvertComposeToBase<true, 10>(inner_function);
#endif
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during ComposeSensitivity. Please help us by "
          "writing an issue about this case at [ "
          "github.com/tataratat/splinepy ]");
    }
  } else {
    switch (inner_function->SplinepyParaDim()) {
    case 1:
      return ConvertComposeToBase<false, 1>(inner_function);
    case 2:
      return ConvertComposeToBase<false, 2>(inner_function);
    case 3:
      return ConvertComposeToBase<false, 3>(inner_function);
#ifdef SPLINEPY_MORE
    case 4:
      return ConvertComposeToBase<false, 4>(inner_function);
    case 5:
      return ConvertComposeToBase<false, 5>(inner_function);
    case 6:
      return ConvertComposeToBase<false, 6>(inner_function);
    case 7:
      return ConvertComposeToBase<false, 7>(inner_function);
    case 8:
      return ConvertComposeToBase<false, 8>(inner_function);
    case 9:
      return ConvertComposeToBase<false, 9>(inner_function);
    case 10:
      return ConvertComposeToBase<false, 10>(inner_function);
#endif
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during ComposeSensitivity. Please help us by "
          "writing an issue about this case at [ "
          "github.com/tataratat/splinepy ]");
    }
  }
  splinepy::utils::PrintAndThrowError(
      "Something went wrong during ComposeSensitivity. Please help us by "
      "writing an issue about this case at [ "
      "github.com/tataratat/splinepy ]");
  // make compiler happy
  return std::vector<std::shared_ptr<SplinepyBase>>{};
}

template<std::size_t para_dim, std::size_t dim>
std::shared_ptr<SplinepyBase> Bezier<para_dim, dim>::SplinepyCompose(
    const std::shared_ptr<SplinepyBase>& inner_function) const {
  // type check
  if (inner_function->SplinepySplineName().find("Bezier")
      == std::string::npos) {
    splinepy::utils::PrintAndThrowError(
        "Bezier composition requires inner function to be a bezier type.",
        "Given inner function -",
        inner_function->SplinepyWhatAmI());
  }

  // composable?
  if (inner_function->SplinepyDim() != this->SplinepyParaDim()) {
    splinepy::utils::PrintAndThrowError(
        "Spline composition requires inner function to have same physical",
        "dimension as outer function's parametric dimension.",
        "Outer Function:",
        this->SplinepyWhatAmI(),
        "/",
        "Inner Function:",
        inner_function->SplinepyWhatAmI());
  }
  // compose - return correct type.
  if (inner_function->SplinepyIsRational()) {
    switch (inner_function->SplinepyParaDim()) {
    case 1:
      return std::make_shared<RationalBezier<1, dim>>(this->Compose(
          static_cast<RationalBezier<1, para_dim>&>(*inner_function)));
    case 2:
      return std::make_shared<RationalBezier<2, dim>>(this->Compose(
          static_cast<RationalBezier<2, para_dim>&>(*inner_function)));
    case 3:
      return std::make_shared<RationalBezier<3, dim>>(this->Compose(
          static_cast<RationalBezier<3, para_dim>&>(*inner_function)));
#ifdef SPLINEPY_MORE
    case 4:
      return std::make_shared<RationalBezier<4, dim>>(this->Compose(
          static_cast<RationalBezier<4, para_dim>&>(*inner_function)));
    case 5:
      return std::make_shared<RationalBezier<5, dim>>(this->Compose(
          static_cast<RationalBezier<5, para_dim>&>(*inner_function)));
    case 6:
      return std::make_shared<RationalBezier<6, dim>>(this->Compose(
          static_cast<RationalBezier<6, para_dim>&>(*inner_function)));
    case 7:
      return std::make_shared<RationalBezier<7, dim>>(this->Compose(
          static_cast<RationalBezier<7, para_dim>&>(*inner_function)));
    case 8:
      return std::make_shared<RationalBezier<8, dim>>(this->Compose(
          static_cast<RationalBezier<8, para_dim>&>(*inner_function)));
    case 9:
      return std::make_shared<RationalBezier<9, dim>>(this->Compose(
          static_cast<RationalBezier<9, para_dim>&>(*inner_function)));
    case 10:
      return std::make_shared<RationalBezier<10, dim>>(this->Compose(
          static_cast<RationalBezier<10, para_dim>&>(*inner_function)));
#endif
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during Compose. Please help us by writing an "
          "issue about this case at [ github.com/tataratat/splinepy ]");
    }
  } else {
    switch (inner_function->SplinepyParaDim()) {
    case 1:
      return std::make_shared<Bezier<1, dim>>(
          this->Compose(static_cast<Bezier<1, para_dim>&>(*inner_function)));
    case 2:
      return std::make_shared<Bezier<2, dim>>(
          this->Compose(static_cast<Bezier<2, para_dim>&>(*inner_function)));
    case 3:
      return std::make_shared<Bezier<3, dim>>(
          this->Compose(static_cast<Bezier<3, para_dim>&>(*inner_function)));
#ifdef SPLINEPY_MORE
    case 4:
      return std::make_shared<Bezier<4, dim>>(
          this->Compose(static_cast<Bezier<4, para_dim>&>(*inner_function)));
    case 5:
      return std::make_shared<Bezier<5, dim>>(
          this->Compose(static_cast<Bezier<5, para_dim>&>(*inner_function)));
    case 6:
      return std::make_shared<Bezier<6, dim>>(
          this->Compose(static_cast<Bezier<6, para_dim>&>(*inner_function)));
    case 7:
      return std::make_shared<Bezier<7, dim>>(
          this->Compose(static_cast<Bezier<7, para_dim>&>(*inner_function)));
    case 8:
      return std::make_shared<Bezier<8, dim>>(
          this->Compose(static_cast<Bezier<8, para_dim>&>(*inner_function)));
    case 9:
      return std::make_shared<Bezier<9, dim>>(
          this->Compose(static_cast<Bezier<9, para_dim>&>(*inner_function)));
    case 10:
      return std::make_shared<Bezier<10, dim>>(
          this->Compose(static_cast<Bezier<10, para_dim>&>(*inner_function)));
#endif
    default:
      splinepy::utils::PrintAndThrowError(
          "Something went wrong during Compose. Please help us by writing an "
          "issue about this case at [ github.com/tataratat/splinepy ]");
    }
  }
  splinepy::utils::PrintAndThrowError(
      "Something went very wrong during Compose. Please help us by writing "
      "an issue about this case at [ github.com/tataratat/splinepy ]");
  // make compiler happy
  return std::shared_ptr<SplinepyBase>{};
}

template<std::size_t para_dim, std::size_t dim>
std::vector<std::shared_ptr<SplinepyBase>>
Bezier<para_dim, dim>::SplinepySplit(const int& p_dim,
                                     const double& location) const {
  // split
  auto bm_splitted =
      Base_::SplitAtPosition(location, static_cast<std::size_t>(p_dim));

  // make it splinepybase
  std::vector<std::shared_ptr<SplinepyBase>> split;
  const std::size_t n_splitted = bm_splitted.size();
  split.reserve(n_splitted); // this should be always 2
  for (std::size_t i{}; i < n_splitted; ++i) {
    split.emplace_back(std::make_shared<Bezier<para_dim, dim>>(bm_splitted[i]));
  }

  return split;
}

template<std::size_t para_dim, std::size_t dim>
std::shared_ptr<SplinepyBase>
Bezier<para_dim, dim>::SplinepyDerivativeSpline(const int* orders) const {
  // copy construct to start
  Base_ derived_bez{*this};
  Derivative_ orders_{};
  std::copy(orders, orders + para_dim, std::begin(orders_));
  // derive
  derived_bez = derived_bez.DerivativeWRTParametricDimension(orders_);

  return std::make_shared<Bezier<para_dim, dim>>(derived_bez);
}

template<std::size_t para_dim, std::size_t dim>
std::vector<std::shared_ptr<SplinepyBase>>
Bezier<para_dim, dim>::SplinepyExtractBezierPatches() const {
  // should copy
  return {std::make_shared<Bezier<para_dim, dim>>(*this)};
}

template<std::size_t para_dim, std::size_t dim>
std::shared_ptr<SplinepyBase>
Bezier<para_dim, dim>::SplinepyExtractBoundary(const int& boundary_id) {
  return splinepy::splines::helpers::ExtractBoundaryMeshSlice(*this,
                                                              boundary_id);
}

template<std::size_t para_dim, std::size_t dim>
std::shared_ptr<SplinepyBase>
Bezier<para_dim, dim>::SplinepyExtractDim(const int& phys_dim) const {
  return std::make_shared<Bezier<para_dim, 1>>(
      Base_::ExtractDimension(static_cast<std::size_t>(phys_dim)));
}

template<std::size_t para_dim, std::size_t dim>
std::shared_ptr<SplinepyBase>
Bezier<para_dim, dim>::SplinepyCompositionDerivative(
    const std::shared_ptr<SplinepyBase>& inner,
    const std::shared_ptr<SplinepyBase>& inner_derivative) const {
  // type check
  if (inner->SplinepySplineName().find("Bezier") == std::string::npos) {
    splinepy::utils::PrintAndThrowError(
        "Bezier composition derivative requires inner function to be",
        "a bezier type. Given inner function -",
        inner->SplinepyWhatAmI());
  }

  // composable?
  if (inner->SplinepyDim() != this->SplinepyParaDim()) {
    splinepy::utils::PrintAndThrowError(
        "Spline composition requires inner function to have same physical",
        "dimension as outer function's parametric dimension.",
        "Outer Function:",
        this->SplinepyWhatAmI(),
        "/",
        "Inner Function:",
        inner->SplinepyWhatAmI());
  }
  std::array<int, para_dim> order_query;
  order_query.fill(0);
  order_query[0] = 1;
  std::shared_ptr<SplinepyBase> composition_derivative =
      this->SplinepyDerivativeSpline(order_query.data())
          ->SplinepyCompose(inner)
          ->SplinepyMultiply(inner_derivative->SplinepyExtractDim(0));

  for (std::size_t i{1}; i < para_dim; ++i) {
    order_query.fill(0);
    order_query[i] = 1;
    composition_derivative = composition_derivative->SplinepyAdd(
        this->SplinepyDerivativeSpline(order_query.data())
            ->SplinepyCompose(inner)
            ->SplinepyMultiply(inner_derivative->SplinepyExtractDim(i)));
  }
  return composition_derivative;
}

} // namespace splinepy::splines

// #include "splinepy/explicit/bezier.hpp"
