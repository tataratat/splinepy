namespace splinepy::splines {

template<std::size_t para_dim>
typename RationalBezier<para_dim>::Base_ RationalBezier<para_dim>::CreateBase(const int* degrees,
                                                  const double* control_points,
                                                  const double* weights,
                                                  const int dim) {

  std::array<std::size_t, para_dim> bm_degrees{};
  std::size_t ncps{1};

  // formulate degrees
  for (std::size_t i{}; i < para_dim; ++i) {
    bm_degrees[i] = degrees[i];
    ncps *= degrees[i] + 1;
  }

  // formulate control_points
  VectorSpace_ bm_control_points{};
  bm_control_points.SetDimension(dim);
  auto& bm_vertices = bm_control_points.GetVertices();
  bm_vertices.reserve(ncps * dim);

  std::vector<Weight_> bm_weights{};
  bm_weights.reserve(ncps);

  for (std::size_t i{}; i < ncps; ++i) {
    const double* this_cp_ptr = &control_points[i * dim];
    for (int j{}; j < dim; ++j) {
      bm_vertices.push_back(this_cp_ptr[j]);
    }
    bm_weights.push_back(weights[i]);
  }

  return Base_(bm_degrees, bm_control_points, bm_weights);
}

/// @brief Function wrapper, also for helper functions. Evaluates query
/// parametric coordinate
/// @param query
template<std::size_t para_dim>
constexpr auto
RationalBezier<para_dim>::operator()(const ParametricCoordinate_& query) const {
  return Base_::Evaluate(query);
}

/// @brief Evaluates derivative of spline
/// @param query
/// @param order
template<std::size_t para_dim>
constexpr auto
RationalBezier<para_dim>::operator()(const ParametricCoordinate_& query,
                                     const Derivative_& order) const {
  return Base_::EvaluateDerivative(query, order);
}

/// @brief Elevate degree
/// @param p_dim
template<std::size_t para_dim>
constexpr auto RationalBezier<para_dim>::ElevateDegree(const Dimension_ p_dim) {
  return Base_::OrderElevateAlongParametricDimension(p_dim);
}

/// @brief Number of supports
/// @return int
template<std::size_t para_dim>
int RationalBezier<para_dim>::SplinepyNumberOfSupports() const {
  return helpers::GetNumberOfSupports(*this);
}

/// @brief Current properties
template<std::size_t para_dim>
void RationalBezier<para_dim>::SplinepyCurrentProperties(
    int* degrees,
    std::vector<std::vector<double>>* knot_vectors /* untouched */,
    double* control_points,
    double* weights) const {

  // degrees
  if (degrees) {
    for (std::size_t i{}; i < kParaDim; ++i) {
      degrees[i] = static_cast<int>(Base_::GetDegrees()[i]);
    }
  }

  // control_points and weights
  if (control_points || weights) {

    const auto& bm_weighted_cps = Base_::GetWeightedControlPoints();
    const auto& bm_weighted_vertices = bm_weighted_cps.GetVertices();
    const auto& bm_weights = Base_::GetWeights(); // this is vector
    const auto ncps = bm_weighted_vertices.size();
    const auto dim = bm_weighted_cps.GetDimension();

    auto vertices_iter = bm_weighted_vertices.cbegin();

    for (std::size_t i{}; i < ncps; ++i) {
      const double& w = bm_weights[i];
      if (weights) {
        weights[i] = w;
      }
      if (control_points) {
        const double inv_weight = static_cast<double>(1.) / w;
        double* this_cp = &control_points[i * dim];
        for (std::size_t j{}; j < dim; ++j) {
          this_cp[j] = *(vertices_iter++) * inv_weight;
        }
      }
    }
  }
}

template<std::size_t para_dim>
void RationalBezier<para_dim>::SplinepyParametricBounds(
    double* para_bounds) const {
  for (std::size_t i{}; i < kParaDim; ++i) {
    // lower bounds
    para_bounds[i] = 0.;
    // upper bounds
    para_bounds[kParaDim + i] = 1.;
  }
}

template<std::size_t para_dim>
void RationalBezier<para_dim>::SplinepyControlMeshResolutions(
    int* control_mesh_res) const {
  const auto cm_res =
      helpers::GetControlMeshResolutions(*this);
  std::copy_n(cm_res.begin(), para_dim, control_mesh_res);
}

/// @brief Calculate Greville abscissae for Rational Bezier
///
/// @param[out] greville_abscissae pointer to solution
/// @param[in] i_para_dim parametric dimension
template<std::size_t para_dim>
void RationalBezier<para_dim>::SplinepyGrevilleAbscissae(
    double* greville_abscissae,
    const int& i_para_dim) const {
  helpers::GetGrevilleAbscissae(*this,
                                                   greville_abscissae,
                                                   i_para_dim);
}

template<std::size_t para_dim>
void RationalBezier<para_dim>::SplinepyEvaluate(const double* para_coord,
                                                double* evaluated) const {
  helpers::ScalarTypeEvaluate(*this, para_coord, evaluated);
}
template<std::size_t para_dim>
void RationalBezier<para_dim>::SplinepyDerivative(const double* para_coord,
                                                  const int* orders,
                                                  double* derived) const {
  helpers::ScalarTypeDerivative(*this,
                                                   para_coord,
                                                   orders,
                                                   derived);
}

template<std::size_t para_dim>
void RationalBezier<para_dim>::SplinepyJacobian(const double* para_coord,
                                                double* jacobians) const {
  helpers::ScalarTypeJacobian(*this, para_coord, jacobians);
}

template<std::size_t para_dim>
void RationalBezier<para_dim>::SplinepyElevateDegree(const int& p_dim) {
  helpers::ScalarTypeElevateDegree(*this, p_dim);
}

template<std::size_t para_dim>
void RationalBezier<para_dim>::SplinepyBasis(const double* para_coord,
                                             double* basis) const {
  helpers::BezierBasis(*this, para_coord, basis);
}

template<std::size_t para_dim>
void RationalBezier<para_dim>::SplinepyBasisDerivative(
    const double* para_coord,
    const int* order,
    double* basis_der) const {
  helpers::BezierBasisDerivative(*this,
                                                    para_coord,
                                                    order,
                                                    basis_der);
}

template<std::size_t para_dim>
void RationalBezier<para_dim>::SplinepySupport(const double* para_coord,
                                               int* support) const {
  helpers::BezierSupport(*this, para_coord, support);
}

/// Basis Function values and their support IDs
template<std::size_t para_dim>
void RationalBezier<para_dim>::SplinepyBasisAndSupport(const double* para_coord,
                                                       double* basis,
                                                       int* support) const {

  SplinepyBasis(para_coord, basis);
  SplinepySupport(para_coord, support);
}

/// Basis Function Derivative and their support IDs
template<std::size_t para_dim>
void RationalBezier<para_dim>::SplinepyBasisDerivativeAndSupport(
    const double* para_coord,
    const int* orders,
    double* basis_der,
    int* support) const {
  SplinepyBasisDerivative(para_coord, orders, basis_der);
  SplinepySupport(para_coord, support);
}

template<std::size_t para_dim>
void RationalBezier<para_dim>::SplinepyPlantNewKdTreeForProximity(
    const int* resolutions,
    const int& nthreads) {
  GetProximity().PlantNewKdTree(resolutions, nthreads);
}

/// Verbose proximity query - make sure to plant a kdtree first.
template<std::size_t para_dim>
void RationalBezier<para_dim>::SplinepyVerboseProximity(
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

/// only applicable to the splines of same para_dim, same type, and
/// {1 or same} dim.
template<std::size_t para_dim>
std::shared_ptr<SplinepyBase> RationalBezier<para_dim>::SplinepyMultiply(
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
    auto true_a = static_cast<RationalBezier<para_dim>&>(*a);
    return std::make_shared<RationalBezier<para_dim>>(
        this->Base_::operator*(true_a));
  } else {
    auto true_a = static_cast<RationalBezier<para_dim>&>(*a);
    return std::make_shared<RationalBezier<para_dim>>(
        this->Base_::operator*(true_a));
  }
}

/// Spline addition.
/// requires same type, para_dim, dim
template<std::size_t para_dim>
std::shared_ptr<SplinepyBase> RationalBezier<para_dim>::SplinepyAdd(
    const std::shared_ptr<SplinepyBase>& a) const {
  SplinepySplineNameMatches(
      *this,
      *a,
      "Spline addition requires splines of the same time.",
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

  return std::make_shared<RationalBezier<para_dim>>(
      this->Base_::operator+(static_cast<RationalBezier<para_dim>&>(*a)));
}

/**
 * @brief Convert the composed splines into base splines for python side
 *
 * @tparam is_rational Flag for rational beziers
 * @tparam inner_para_dim parametric of inner function
 * @param inner_function inner function of composition
 * @return std::vector<std::shared_ptr<SplinepyBase>> Vector of SplinepyBase
 * splines
 */
template<std::size_t para_dim>
template<bool is_rational, size_t inner_para_dim>
std::vector<std::shared_ptr<SplinepyBase>>
RationalBezier<para_dim>::ConvertComposeToBase(
    const std::shared_ptr<SplinepyBase>& inner_function) const {
  if constexpr (is_rational) {
    auto super_return = this->ComposeSensitivity(
        static_cast<RationalBezier<inner_para_dim>&>(*inner_function));

    std::vector<std::shared_ptr<SplinepyBase>> to_return(super_return.size());

    for (size_t i{}; i < super_return.size(); ++i) {
      to_return[i] =
          std::make_shared<RationalBezier<inner_para_dim>>(super_return[i]);
    }
    return to_return;
  } else {
    auto super_return = this->ComposeSensitivity(
        static_cast<Bezier<inner_para_dim>&>(*inner_function));

    std::vector<std::shared_ptr<SplinepyBase>> to_return(super_return.size());

    for (size_t i{}; i < super_return.size(); ++i) {
      to_return[i] =
          std::make_shared<RationalBezier<inner_para_dim>>(super_return[i]);
    }
    return to_return;
  }
}

/**
 * @brief Compute sensitivities of Composed splines with respect to outer
 * functions control points
 *
 *  Spline composition.
 *  inner_function requirements:
 *    1. Bezier Types
 *    2. dim is same as outer_function's par_dim
 *
 * @param inner_function inner function of composition
 * @return std::vector<std::shared_ptr<SplinepyBase>>
 */
template<std::size_t para_dim>
std::vector<std::shared_ptr<SplinepyBase>>
RationalBezier<para_dim>::SplinepyComposeSensitivities(
    const std::shared_ptr<SplinepyBase>& inner_function) const {
  // type check
  if (inner_function->SplinepySplineName().find("Bezier")
      == std::string::npos) {
    utils::PrintAndThrowError(
        "Bezier composition requires inner function to be a bezier type.",
        "Given inner function -",
        inner_function->SplinepyWhatAmI());
  }

  // composable?
  if (inner_function->SplinepyDim() != this->SplinepyParaDim()) {
    utils::PrintAndThrowError(
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
      utils::PrintAndThrowError(
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
      utils::PrintAndThrowError(
          "Something went wrong during ComposeSensitivity. Please help us by "
          "writing an issue about this case at [ "
          "github.com/tataratat/splinepy ]");
    }
  }
  utils::PrintAndThrowError(
      "Something went wrong during ComposeSensitivity. Please help us by "
      "writing an issue about this case at [ "
      "github.com/tataratat/splinepy ]");
  // make compiler happy
  return std::vector<std::shared_ptr<SplinepyBase>>{};
}

/// Spline composition.
/// inner_function requirements:
///   1. Bezier Types
///   2. dim is same as outer_function's par_dim
template<std::size_t para_dim>
std::shared_ptr<SplinepyBase> RationalBezier<para_dim>::SplinepyCompose(
    const std::shared_ptr<SplinepyBase>& inner_function) const {
  // type check
  if (inner_function->SplinepySplineName().find("Bezier")
      == std::string::npos) {
    utils::PrintAndThrowError(
        "Bezier composition requires inner function to be a bezier type.",
        "Given inner function -",
        inner_function->SplinepyWhatAmI());
  }

  // composable?
  if (inner_function->SplinepyDim() != this->SplinepyParaDim()) {
    utils::PrintAndThrowError(
        "Spline composition requires inner function to have same physical",
        "dimension as outer function's parametric dimension.",
        "Outer Function:",
        this->SplinepyWhatAmI(),
        "/",
        "Inner Function:",
        inner_function->SplinepyWhatAmI());
  }

  // compose - this time, it is always rational.
  if (inner_function->SplinepyIsRational()) {
    switch (inner_function->SplinepyParaDim()) {
    case 1:
      return std::make_shared<RationalBezier<1>>(
          this->Compose(static_cast<RationalBezier<1>&>(*inner_function)));
    case 2:
      return std::make_shared<RationalBezier<2>>(
          this->Compose(static_cast<RationalBezier<2>&>(*inner_function)));
    case 3:
      return std::make_shared<RationalBezier<3>>(
          this->Compose(static_cast<RationalBezier<3>&>(*inner_function)));
#ifdef SPLINEPY_MORE
    case 4:
      return std::make_shared<RationalBezier<4>>(
          this->Compose(static_cast<RationalBezier<4>&>(*inner_function)));
    case 5:
      return std::make_shared<RationalBezier<5>>(
          this->Compose(static_cast<RationalBezier<5>&>(*inner_function)));
    case 6:
      return std::make_shared<RationalBezier<6>>(
          this->Compose(static_cast<RationalBezier<6>&>(*inner_function)));
    case 7:
      return std::make_shared<RationalBezier<7>>(
          this->Compose(static_cast<RationalBezier<7>&>(*inner_function)));
    case 8:
      return std::make_shared<RationalBezier<8>>(
          this->Compose(static_cast<RationalBezier<8>&>(*inner_function)));
    case 9:
      return std::make_shared<RationalBezier<9>>(
          this->Compose(static_cast<RationalBezier<9>&>(*inner_function)));
    case 10:
      return std::make_shared<RationalBezier<10>>(
          this->Compose(static_cast<RationalBezier<10>&>(*inner_function)));
#endif
    default:
      utils::PrintAndThrowError(
          "Something went wrong during Compose. Please help us by writing an "
          "issue about this case at [ github.com/tataratat/splinepy ]");
    }
  } else {
    switch (inner_function->SplinepyParaDim()) {
    case 1:
      return std::make_shared<RationalBezier<1>>(
          this->Compose(static_cast<Bezier<1>&>(*inner_function)));
    case 2:
      return std::make_shared<RationalBezier<2>>(
          this->Compose(static_cast<Bezier<2>&>(*inner_function)));
    case 3:
      return std::make_shared<RationalBezier<3>>(
          this->Compose(static_cast<Bezier<3>&>(*inner_function)));
#ifdef SPLINEPY_MORE
    case 4:
      return std::make_shared<RationalBezier<4>>(
          this->Compose(static_cast<Bezier<4>&>(*inner_function)));
    case 5:
      return std::make_shared<RationalBezier<5>>(
          this->Compose(static_cast<Bezier<5>&>(*inner_function)));
    case 6:
      return std::make_shared<RationalBezier<6>>(
          this->Compose(static_cast<Bezier<6>&>(*inner_function)));
    case 7:
      return std::make_shared<RationalBezier<7>>(
          this->Compose(static_cast<Bezier<7>&>(*inner_function)));
    case 8:
      return std::make_shared<RationalBezier<8>>(
          this->Compose(static_cast<Bezier<8>&>(*inner_function)));
    case 9:
      return std::make_shared<RationalBezier<9>>(
          this->Compose(static_cast<Bezier<9>&>(*inner_function)));
    case 10:
      return std::make_shared<RationalBezier<10>>(
          this->Compose(static_cast<Bezier<10>&>(*inner_function)));
#endif
    default:
      utils::PrintAndThrowError(
          "Something went wrong during Compose. Please help us by writing an "
          "issue about this case at [ github.com/tataratat/splinepy ]");
    }
  }
  utils::PrintAndThrowError(
      "Something went very wrong during Compose. Please help us by writing"
      "an issue about this case at [ github.com/tataratat/splinepy ]");
  // make compiler happy
  return std::shared_ptr<SplinepyBase>{};
}

template<std::size_t para_dim>
std::vector<std::shared_ptr<SplinepyBase>>
RationalBezier<para_dim>::SplinepySplit(const int& p_dim,
                                        const double& location) const {
  // split
  auto bm_splitted =
      Base_::SplitAtPosition(location, static_cast<std::size_t>(p_dim));

  // make it splinepybase
  std::vector<std::shared_ptr<SplinepyBase>> splitted;
  const std::size_t n_splitted = bm_splitted.size();
  splitted.reserve(n_splitted); // this should be always 2
  for (std::size_t i{}; i < n_splitted; ++i) {
    splitted.emplace_back(
        std::make_shared<RationalBezier<para_dim>>(bm_splitted[i]));
  }

  return splitted;
}

template<std::size_t para_dim>
std::shared_ptr<SplinepyBase>
RationalBezier<para_dim>::SplinepyDerivativeSpline(const int* orders) const {
  // copy construct to start
  Base_ derived_bez{*this};
  // derive
  for (std::size_t i{}; i < para_dim; ++i) {
    for (int j{}; j < orders[i]; ++j) {
      derived_bez = derived_bez.DerivativeWRTParametricDimension(i);
    }
  }

  return std::make_shared<RationalBezier<para_dim>>(derived_bez);
}

template<std::size_t para_dim>
std::shared_ptr<SplinepyBase>
RationalBezier<para_dim>::SplinepyExtractBoundary(const int& boundary_id) {
  return helpers::ExtractBoundaryMeshSlice(*this,
                                                              boundary_id);
}

template<std::size_t para_dim>
std::vector<std::shared_ptr<SplinepyBase>>
RationalBezier<para_dim>::SplinepyExtractBezierPatches() const {
  // should copy
  return {std::make_shared<RationalBezier<para_dim>>(*this)};
}

template<std::size_t para_dim>
std::shared_ptr<SplinepyBase>
RationalBezier<para_dim>::SplinepyExtractDim(const int& phys_dim) const {
  return std::make_shared<RationalBezier<para_dim>>(
      Base_::ExtractDimension(static_cast<std::size_t>(phys_dim)));
}

/// Derivative of composition
template<std::size_t para_dim>
std::shared_ptr<SplinepyBase>
RationalBezier<para_dim>::SplinepyCompositionDerivative(
    const std::shared_ptr<SplinepyBase>& inner,
    const std::shared_ptr<SplinepyBase>& inner_derivative) const {
  // type check
  if (inner->SplinepySplineName().find("Bezier") == std::string::npos) {
    utils::PrintAndThrowError(
        "Bezier composition derivative requires inner function to be",
        "a bezier type. Given inner function -",
        inner->SplinepyWhatAmI());
  }

  // composable?
  if (inner->SplinepyDim() != this->SplinepyParaDim()) {
    utils::PrintAndThrowError(
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

template<std::size_t para_dim>
constexpr void RationalBezier<para_dim>::CreateProximityIfNull() {
  if (!proximity_) {
    switch (SplinepyDim()) {
    case 1:
      proximity_ = std::make_unique<Proximity_<1>>(*this);
      break;
    case 2:
      proximity_ = std::make_unique<Proximity_<2>>(*this);
      break;
    case 3:
      proximity_ = std::make_unique<Proximity_<3>>(*this);
      break;
    case 4:
      proximity_ = std::make_unique<Proximity_<4>>(*this);
      break;
    case 5:
      proximity_ = std::make_unique<Proximity_<5>>(*this);
      break;
    case 6:
      proximity_ = std::make_unique<Proximity_<6>>(*this);
      break;
    case 7:
      proximity_ = std::make_unique<Proximity_<7>>(*this);
      break;
    case 8:
      proximity_ = std::make_unique<Proximity_<8>>(*this);
      break;
    case 9:
      proximity_ = std::make_unique<Proximity_<9>>(*this);
      break;
    case 10:
      proximity_ = std::make_unique<Proximity_<10>>(*this);
      break;
    default:
      utils::PrintAndThrowError(
          "Something went wrong during Proximity. Please help us by writing "
          "an "
          "issue about this case at [ github.com/tataratat/splinepy ]");
    }
  }
}

/// @brief Gets proximity
template<std::size_t para_dim>
constexpr typename RationalBezier<para_dim>::ProximityBase_& RationalBezier<para_dim>::GetProximity() {
  CreateProximityIfNull();
  return *proximity_;
}
/// @brief Gets proximity
template<std::size_t para_dim>
constexpr const typename RationalBezier<para_dim>::ProximityBase_& RationalBezier<para_dim>::GetProximity() const {
  return *proximity_;
}

}