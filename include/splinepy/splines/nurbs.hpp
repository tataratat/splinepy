#pragma once

// SplineLib
#include <BSplineLib/Splines/nurbs.hpp>

#include <splinepy/proximity/proximity.hpp>
#include <splinepy/splines/helpers/basis_functions.hpp>
#include <splinepy/splines/helpers/extract.hpp>
#include <splinepy/splines/helpers/extract_bezier_patches.hpp>
#include <splinepy/splines/helpers/properties.hpp>
#include <splinepy/splines/helpers/scalar_type_wrapper.hpp>
#include <splinepy/splines/splinepy_base.hpp>

namespace splinepy::splines {

/// @class Nurbs
/// @brief Non-uniform rational B-spline (NURBS) class
/// @tparam para_dim Dimension of parametric space
/// @tparam dim Dimension of physical space
template<int para_dim, int dim>
class Nurbs : public splinepy::splines::SplinepyBase,
              public bsplinelib::splines::Nurbs<para_dim, dim> {
public:
  /// @brief Dimension of parametric space
  static constexpr int kParaDim = para_dim;
  /// @brief Dimension of physical space
  static constexpr int kDim = dim;
  /// @brief It is a rational spline
  static constexpr bool kIsRational = true;
  /// @brief It has knot vectors
  static constexpr bool kHasKnotVectors = true;

  // TODO remove after test
  /// @brief Dimension of parameter space
  constexpr static int para_dim_ = para_dim;
  /// @brief Dimension of physical space
  constexpr static int dim_ = dim;

  // self
  template<int s_para_dim, int s_dim>
  using SelfTemplate_ = Nurbs<s_para_dim, s_dim>;

  // splinepy
  using SplinepyBase_ = splinepy::splines::SplinepyBase;
  using WeightedControlPointPointers_ =
      typename SplinepyBase_::WeightedControlPointPointers_;
  using WeightPointers_ = typename SplinepyBase_::WeightPointers_;

  // splinelib
  using Base_ = bsplinelib::splines::Nurbs<para_dim, dim>;
  // Parameter space
  using ParameterSpace_ = typename Base_::ParameterSpace_;
  using Degrees_ = typename ParameterSpace_::Degrees_;
  using Degree_ = typename Degrees_::value_type;
  using KnotVectors_ = typename ParameterSpace_::KnotVectors_;
  using KnotVector_ = typename KnotVectors_::value_type::element_type;
  using Knots_ = typename Base_::Base_::Knots_;
  using Knot_ = typename Base_::Knot_;
  using KnotRatios_ = typename Base_::ParameterSpace_::KnotRatios_;
  using KnotRatio_ = typename KnotRatios_::value_type;
  using ParametricCoordinate_ = typename Base_::ParametricCoordinate_;
  using ScalarParametricCoordinate_ =
      typename ParametricCoordinate_::value_type;
  // Weighted Vector Space
  using WeightedVectorSpace_ = typename Base_::WeightedVectorSpace_;
  using PhysicalSpace_ = WeightedVectorSpace_;
  using Coordinates_ = typename WeightedVectorSpace_::Coordinates_;
  using HomogeneousCoordinates_ =
      typename WeightedVectorSpace_::Base_::Coordinates_;
  using Coordinate_ = typename Base_::Coordinate_;
  using ScalarCoordinate_ = typename Coordinate_::value_type;
  using Weights_ = typename WeightedVectorSpace_::Weights_;
  using Weight_ = typename Weights_::value_type;
  //  Frequently used types
  using Derivative_ = typename Base_::Derivative_;
  using Dimension_ = bsplinelib::Dimension;
  using Tolerance_ = bsplinelib::splines::Tolerance;
  using OutputInformation_ =
      bsplinelib::Tuple<typename ParameterSpace_::OutputInformation_,
                        typename WeightedVectorSpace_::OutputInformation_>;
  using Index_ = typename Base_::Base_::Index_;
  using IndexLength_ = typename Index_::Length_;
  using IndexValue_ = typename Index_::Value_;
  // Advanced use
  using HomogeneousBSpline_ = typename Base_::HomogeneousBSpline_;
  using Proximity_ = splinepy::proximity::Proximity<Nurbs<para_dim, dim>>;

  /// @brief raw ptr based inithelper.
  /// @param degrees should have same size as parametric dimension
  /// @param knot_vectors as vector of vector, we can keep track of their
  /// length, as well as the legnth of control_points/weights.
  /// @param control_points
  /// @param weights
  static Base_ CreateBase(const int* degrees,
                          const std::vector<std::vector<double>>& knot_vectors,
                          const double* control_points,
                          const double* weights) {
    // process all the info and turn them into SplineLib types to initialize
    // Base_.

    // Prepare temporary containers
    Degrees_ sl_degrees;            // std::array
    Knots_ sl_knots;                // std::vector
    KnotVectors_ sl_knot_vectors;   // std::array
    Coordinates_ sl_control_points; // std::vector
    Weights_ sl_weights;            // std::vector
    std::size_t ncps{1};

    // Formulate degrees and knotvectors
    for (std::size_t i{}; i < kParaDim; ++i) {
      // degrees
      sl_degrees[i] = Degree_(degrees[i]);
      // knot vectors
      const auto& knot_vector = knot_vectors[i];
      std::size_t nkv = knot_vector.size();
      sl_knots.clear();
      sl_knots.reserve(nkv);
      // try if this works after namedtype ext
      // sl_knots = knot_vector;
      for (std::size_t j{}; j < nkv; ++j) {
        sl_knots.emplace_back(Knot_{knot_vector[j]});
      }
      std::shared_ptr sl_knot_vector{std::make_shared<KnotVector_>(sl_knots)};
      sl_knot_vectors[i] = sl_knot_vector;

      ncps *= nkv - degrees[i] - 1;
    }

    // Formulate ParameterSpace
    auto sl_parameter_space =
        std::make_shared<ParameterSpace_>(sl_knot_vectors, sl_degrees);

    // Formulate control_points and weights
    sl_control_points.reserve(ncps);
    sl_weights.reserve(ncps);
    for (std::size_t i{}; i < ncps; ++i) {
      // control_points
      Coordinate_ control_point;
      for (std::size_t j{}; j < kDim; ++j) {
        control_point[j] = ScalarCoordinate_{control_points[i * kDim + j]};
      }
      sl_control_points.push_back(std::move(control_point));

      // weights
      sl_weights.emplace_back(weights[i]);
    }

    auto sl_weighted_space =
        std::make_shared<WeightedVectorSpace_>(sl_control_points, sl_weights);

    // return init
    return Base_(sl_parameter_space, sl_weighted_space);
  }

  /// @brief Constructor based on raw pointer
  /// @param degrees
  /// @param knot_vectors
  /// @param control_points
  /// @param weights
  Nurbs(const int* degrees,
        const std::vector<std::vector<double>>& knot_vectors,
        const double* control_points,
        const double* weights)
      : Base_(CreateBase(degrees, knot_vectors, control_points, weights)) {}
  // inherit ctor
  using Base_::Base_;

  /// @brief Get degrees of parameter space
  constexpr const Degrees_& GetDegrees() const {
    return GetParameterSpace().GetDegrees();
  };

  /// @brief Get knot vectors
  constexpr const KnotVectors_& GetKnotVectors() const {
    return GetParameterSpace().GetKnotVectors();
  }

  /// @brief Get coordinates. Nurbs' coordinates are homogeneous
  constexpr const HomogeneousCoordinates_& GetCoordinates() const {
    return GetWeightedVectorSpace().GetCoordinates();
  }

  // required implementations
  /// @copydoc splinepy::splines::SplinepyBase::SplinepyParaDim
  virtual int SplinepyParaDim() const { return kParaDim; }

  /// @copydoc splinepy::splines::SplinepyBase::SplinepyDim
  virtual int SplinepyDim() const { return kDim; }

  /// @copydoc splinepy::splines::SplinepyBase::SplinepySplineName
  virtual std::string SplinepySplineName() const { return "NURBS"; }

  /// @copydoc splinepy::splines::SplinepyBase::SplinepyWhatAmI
  virtual std::string SplinepyWhatAmI() const {
    return "NURBS, parametric dimension: " + std::to_string(SplinepyParaDim())
           + ", physical dimension: " + std::to_string(SplinepyDim());
  }

  /// @copydoc splinepy::splines::SplinepyBase::SplinepyHasKnotVectors
  virtual bool SplinepyHasKnotVectors() const { return kHasKnotVectors; }

  /// @copydoc splinepy::splines::SplinepyBase::SplinepyIsRational
  virtual bool SplinepyIsRational() const { return kIsRational; }

  /// @copydoc splinepy::splines::SplinepyBase::SplinepyNumberOfControlPoints
  virtual int SplinepyNumberOfControlPoints() const {
    return GetWeightedVectorSpace().GetNumberOfCoordinates();
  }

  /// @copydoc splinepy::splines::SplinepyBase::SplinepyNumberOfSupports
  virtual int SplinepyNumberOfSupports() const {
    return splinepy::splines::helpers::GetNumberOfSupports(*this);
  }

  /// @copydoc splinepy::splines::SplinepyBase::SplinepyCurrentProperties
  virtual void
  SplinepyCurrentProperties(int* degrees,
                            std::vector<std::vector<double>>* knot_vectors,
                            double* control_points,
                            double* weights) const {
    const auto& parameter_space = GetParameterSpace();
    const auto& vector_space = GetWeightedVectorSpace();

    // degrees
    if (degrees) {
      for (std::size_t i{}; i < kParaDim; ++i) {
        degrees[i] = static_cast<int>(parameter_space.GetDegrees()[i]);
      }
    }

    // knot_vectors
    if (knot_vectors) {
      const auto& core_kvs = parameter_space.GetKnotVectors();
      knot_vectors->clear();
      knot_vectors->reserve(kParaDim);
      for (std::size_t i{}; i < kParaDim; ++i) {
        const auto& core_kv = *core_kvs[i];
        const std::size_t kvsize = static_cast<std::size_t>(core_kv.GetSize());
        std::vector<double> kv;
        kv.reserve(kvsize);
        // Index wants int
        for (int j{}; j < static_cast<int>(kvsize); ++j) {
          kv.emplace_back(static_cast<double>(core_kv[bsplinelib::Index{j}]));
        }
        knot_vectors->push_back(std::move(kv));
      }
    }

    // control_points and weights
    if (control_points || weights) {
      std::size_t ncps = vector_space.GetNumberOfCoordinates();
      // fill it up, phil!
      for (std::size_t i{}; i < ncps; ++i) {
        auto const& coord_named_phil = vector_space[bsplinelib::Index{
            static_cast<bsplinelib::Index::Type_>(i)}];
        // unweight - phil needs to first project before filling.
        if (control_points) {
          auto const& projected_phil =
              WeightedVectorSpace_::Project(coord_named_phil);
          for (std::size_t j{}; j < kDim; ++j) {
            control_points[i * kDim + j] =
                static_cast<double>(projected_phil[j]);
          }
        }
        if (weights) {
          weights[i] = static_cast<double>(coord_named_phil[dim]);
        }
      }
    }
  }

  virtual std::shared_ptr<bsplinelib::parameter_spaces::KnotVector>
  SplinepyKnotVector(const int p_dim) {
    if (!(p_dim < para_dim)) {
      splinepy::utils::PrintAndThrowError(
          "Invalid parametric dimension. Should be smaller than",
          para_dim);
    }
    return Base_::Base_::parameter_space_->GetKnotVectors()[p_dim];
  };

  virtual std::shared_ptr<WeightedControlPointPointers_>
  SplinepyWeightedControlPointPointers() {
    if (SplinepyBase_::control_point_pointers_) {
      return SplinepyBase_::control_point_pointers_;
    }
    // create weighted cps
    auto wcpp = std::make_shared<WeightedControlPointPointers_>();
    wcpp->dim_ = kDim;
    wcpp->for_rational_ = kIsRational;
    wcpp->coordinate_begins_.reserve(SplinepyNumberOfControlPoints());
    // create weights
    auto w = std::make_shared<WeightPointers_>();
    w->weights_.reserve(SplinepyNumberOfControlPoints());

    for (auto& w_control_point :
         Base_::weighted_vector_space_->GetCoordinates()) {
      auto* coord_begin = w_control_point.data();
      wcpp->coordinate_begins_.push_back(coord_begin);
      w->weights_.push_back(coord_begin + kDim);
    }

    // reference each other
    w->control_point_pointers_ = wcpp; // weak-ref
    wcpp->weight_pointers_ = w;

    // save
    SplinepyBase_::control_point_pointers_ = wcpp;

    return wcpp;
  }

  virtual std::shared_ptr<WeightPointers_> SplinepyWeightPointers() {
    return SplinepyWeightedControlPointPointers()->weight_pointers_;
  }

  /// @copydoc splinepy::splines::SplinepyBase::SplinepyParametricBounds
  virtual void SplinepyParametricBounds(double* para_bounds) const {
    const auto pbounds = splinepy::splines::helpers::GetParametricBounds(*this);
    for (std::size_t i{}; i < kParaDim; ++i) {
      // lower bounds
      para_bounds[i] = pbounds[0][i];
      // upper bounds
      para_bounds[kParaDim + i] = pbounds[1][i];
    }
  }

  virtual void SplinepyControlMeshResolutions(int* control_mesh_res) const {
    const auto cm_res =
        splinepy::splines::helpers::GetControlMeshResolutions(*this);
    std::copy_n(cm_res.begin(), para_dim, control_mesh_res);
  }

  /// @brief Calculate Greville abscissae for Nurbs
  ///
  /// @param[out] greville_abscissae pointer to solution
  /// @param[in] i_para_dim parametric dimension
  virtual void SplinepyGrevilleAbscissae(double* greville_abscissae,
                                         const int& i_para_dim) const {
    splinepy::splines::helpers::GetGrevilleAbscissae(*this,
                                                     greville_abscissae,
                                                     i_para_dim);
  }

  virtual void SplinepyEvaluate(const double* para_coord,
                                double* evaluated) const {
    splinepy::splines::helpers::ScalarTypeEvaluate(*this,
                                                   para_coord,
                                                   evaluated);
  }

  virtual void SplinepyDerivative(const double* para_coord,
                                  const int* orders,
                                  double* derived) const {
    splinepy::splines::helpers::ScalarTypeDerivative(*this,
                                                     para_coord,
                                                     orders,
                                                     derived);
  }

  virtual void SplinepyJacobian(const double* para_coord,
                                double* jacobians) const {
    splinepy::splines::helpers::ScalarTypeJacobian(*this,
                                                   para_coord,
                                                   jacobians);
  }

  virtual void SplinepyBasis(const double* para_coord, double* basis) const {
    splinepy::splines::helpers::BSplineBasis(*this, para_coord, basis);
  }

  virtual void SplinepyBasisDerivative(const double* para_coord,
                                       const int* order,
                                       double* basis_der) const {
    splinepy::splines::helpers::BSplineBasisDerivative(*this,
                                                       para_coord,
                                                       order,
                                                       basis_der);
  }

  virtual void SplinepySupport(const double* para_coord, int* support) const {
    splinepy::splines::helpers::BSplineSupport(*this, para_coord, support);
  }

  /// Basis Function values and their support IDs
  virtual void SplinepyBasisAndSupport(const double* para_coord,
                                       double* basis,
                                       int* support) const {

    SplinepyBasis(para_coord, basis);
    SplinepySupport(para_coord, support);
  }

  /// Basis Function Derivative and their support IDs
  virtual void SplinepyBasisDerivativeAndSupport(const double* para_coord,
                                                 const int* orders,
                                                 double* basis_der,
                                                 int* support) const {
    SplinepyBasisDerivative(para_coord, orders, basis_der);
    SplinepySupport(para_coord, support);
  }

  virtual void SplinepyPlantNewKdTreeForProximity(const int* resolutions,
                                                  const int& nthreads) {
    splinepy::splines::helpers::ScalarTypePlantNewKdTreeForProximity(
        *this,
        resolutions,
        nthreads);
  }

  /// Verbose proximity query - make sure to plant a kdtree first.
  virtual void SplinepyVerboseProximity(const double* query,
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

  virtual void SplinepyElevateDegree(const int& p_dim) {
    splinepy::splines::helpers::ScalarTypeElevateDegree(*this, p_dim);
  }

  virtual bool SplinepyReduceDegree(const int& p_dim, const double& tolerance) {
    return splinepy::splines::helpers::ScalarTypeReduceDegree(*this,
                                                              p_dim,
                                                              tolerance);
  }

  virtual bool SplinepyInsertKnot(const int& p_dim, const double& knot) {
    return splinepy::splines::helpers::ScalarTypeInsertKnot(*this, p_dim, knot);
  }

  virtual bool SplinepyRemoveKnot(const int& p_dim,
                                  const double& knot,
                                  const double& tolerance) {
    return splinepy::splines::helpers::ScalarTypeRemoveKnot(*this,
                                                            p_dim,
                                                            knot,
                                                            tolerance);
  }

  virtual std::shared_ptr<SplinepyBase>
  SplinepyExtractBoundary(const int& boundary_id) {
    return splinepy::splines::helpers::ExtractBoundaryMeshSlice(*this,
                                                                boundary_id);
  }

  /// Bezier patch extraction
  virtual std::vector<std::shared_ptr<SplinepyBase>>
  SplinepyExtractBezierPatches() const {
    return splinepy::splines::helpers::ExtractBezierPatches<true>(*this);
  }

  /// @brief Gets parameter space
  constexpr const ParameterSpace_& GetParameterSpace() const {
    return *Base_::Base_::parameter_space_;
  }

  /// @brief Gets weighted vector space
  constexpr const WeightedVectorSpace_& GetWeightedVectorSpace() const {
    return *Base_::weighted_vector_space_;
  }

  /// @brief Gets proximity
  constexpr Proximity_& GetProximity() { return *proximity_; }
  /// @brief Gets proximity
  constexpr const Proximity_& GetProximity() const { return *proximity_; }

protected:
  /// @brief Unique pointer to proximity
  std::unique_ptr<Proximity_> proximity_ = std::make_unique<Proximity_>(*this);

}; /* class Nurbs */

} /* namespace splinepy::splines */

#include "splinepy/explicit/nurbs.hpp"
