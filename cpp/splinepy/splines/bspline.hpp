#pragma once

// SplineLib
#include <Sources/Splines/b_spline.hpp>

#include <splinepy/explicit/splinelib/b_spline_extern.hpp>
#include <splinepy/proximity/proximity.hpp>
#include <splinepy/splines/helpers/extract.hpp>
#include <splinepy/splines/helpers/properties.hpp>
#include <splinepy/splines/helpers/scalar_type_wrapper.hpp>
#include <splinepy/splines/splinepy_base.hpp>

namespace splinepy::splines {

template<int para_dim, int dim>
class BSpline : public splinepy::splines::SplinepyBase,
                public splinelib::sources::splines::BSpline<para_dim, dim> {
public:
  static constexpr int kParaDim = para_dim;
  static constexpr int kDim = dim;
  static constexpr bool kIsRational = false;
  static constexpr bool kHasKnotVectors = true;

  // TODO remve afterwards
  constexpr static int para_dim_ = para_dim;
  constexpr static int dim_ = dim;

  // splinepy
  using SplinepyBase_ = splinepy::splines::SplinepyBase;

  // splinelib
  using Base_ = splinelib::sources::splines::BSpline<para_dim, dim>;
  // parameter space
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
  // vector space
  using VectorSpace_ = typename Base_::VectorSpace_;
  using Coordinates_ = typename VectorSpace_::Coordinates_;
  using Coordinate_ = typename Base_::Coordinate_;
  using ScalarCoordinate_ = typename Coordinate_::value_type;
  // frequently used types
  using Derivative_ = typename Base_::Derivative_;
  using Dimension_ = splinelib::Dimension;
  using Tolerance_ = splinelib::sources::splines::Tolerance;
  using OutputInformation_ =
      splinelib::Tuple<typename ParameterSpace_::OutputInformation_,
                       typename VectorSpace_::OutputInformation_>;
  using Index_ = typename Base_::Base_::Index_;
  using IndexLength_ = typename Index_::Length_;
  using IndexValue_ = typename Index_::Value_;
  // Some private ones. Here we make it public
  using BezierInformation_ =
      typename Base_::ParameterSpace_::BezierInformation_;
  using BinomialRatios_ = typename Base_::ParameterSpace_::BinomialRatios_;
  using BinomialRatio_ = typename BinomialRatios_::value_type;
  // Advanced use
  using Proximity_ = splinepy::proximity::Proximity<BSpline<para_dim, dim>>;

  /** raw ptr based inithelper.
   *  degrees should have same size as parametric dimension
   *  having knot_vectors vector of vector, we can keep track of their length,
   * as well as the legnth of control_points/weights.
   */
  static Base_ CreateBase(const int* degrees,
                          const std::vector<std::vector<double>> knot_vectors,
                          const double* control_points) {
    // process all the info and turn them into SplineLib types to initialize
    // Base_.

    // Prepare temporary containers
    Degrees_ sl_degrees;            // std::array
    Knots_ sl_knots;                // std::vector
    KnotVectors_ sl_knot_vectors;   // std::array
    Coordinates_ sl_control_points; // std::vector
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
    for (std::size_t i{}; i < ncps; ++i) {
      // control_points
      Coordinate_ control_point;
      for (std::size_t j{}; j < kDim; ++j) {
        control_point[j] = ScalarCoordinate_{control_points[i * kDim + j]};
      }
      sl_control_points.push_back(std::move(control_point));
    }

    auto sl_vector_space = std::make_shared<VectorSpace_>(sl_control_points);

    // return init
    return Base_(sl_parameter_space, sl_vector_space);
  }

  // rawptr based ctor
  BSpline(const int* degrees,
          const std::vector<std::vector<double>>& knot_vectors,
          const double* control_points)
      : Base_(CreateBase(degrees, knot_vectors, control_points)) {}
  // inherited ctor
  using Base_::Base_;

  constexpr const Degrees_& GetDegrees() const {
    return GetParameterSpace().GetDegrees();
  };

  // required implementations
  virtual int SplinepyParaDim() const { return kParaDim; }

  virtual int SplinepyDim() const { return kDim; }

  virtual std::string SplinepySplineName() const { return "BSpline"; }

  virtual std::string SplinepyWhatAmI() const {
    return "BSpline, parametric dimension: " + std::to_string(SplinepyParaDim())
           + ", physical dimension: " + std::to_string(SplinepyDim());
  }

  virtual bool SplinepyHasKnotVectors() const { return kHasKnotVectors; }

  virtual bool SplinepyIsRational() const { return kIsRational; }

  virtual int SplinepyNumberOfControlPoints() const {
    return GetVectorSpace().GetNumberOfCoordinates();
  }

  virtual int SplinepyNumberOfSupports() const {
    return splinepy::splines::helpers::GetNumberOfSupports(*this);
  }

  virtual void
  SplinepyCurrentProperties(int* degrees,
                            std::vector<std::vector<double>>* knot_vectors,
                            double* control_points,
                            double* weights /* untouched */) const {

    const auto& parameter_space = GetParameterSpace();
    const auto& vector_space = GetVectorSpace();

    // degrees
    for (std::size_t i{}; i < kParaDim; ++i) {
      degrees[i] = static_cast<int>(parameter_space.GetDegrees()[i]);
    }

    // knot_vectors
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
        kv.emplace_back(static_cast<double>(core_kv[splinelib::Index{j}]));
      }
      knot_vectors->push_back(std::move(kv));
    }

    // control_points and weights
    std::size_t ncps = vector_space.GetNumberOfCoordinates();
    // fill it up, phil!
    for (std::size_t i{}; i < ncps; ++i) {
      auto const& coord_named_phil = vector_space[splinelib::Index{
          static_cast<splinelib::Index::Type_>(i)}];
      for (std::size_t j{}; j < kDim; ++j) {
        control_points[i * kDim + j] = static_cast<double>(coord_named_phil[j]);
      }
    }
  }

  virtual void SplinepyParametricBounds(double* para_bounds) const {
    const auto pbounds = splinepy::splines::helpers::GetParametricBounds(*this);
    for (std::size_t i{}; i < kParaDim; ++i) {
      // lower bounds
      para_bounds[i] = pbounds[0][i];
      // upper bounds
      para_bounds[kParaDim + i] = pbounds[1][i];
    }
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

  virtual void SplinepyBasisAndSupport(const double* para_coord,
                                       double* basis,
                                       int* support) const {
    ParameterSpace_ const& parameter_space = *Base_::Base_::parameter_space_;

    typename ParameterSpace_::UniqueEvaluations_ unique_evaluations;
    parameter_space.template InitializeUniqueEvaluations<false>(
        unique_evaluations);

    ParametricCoordinate_ sl_para_coord;
    for (std::size_t i{}; i < kParaDim; ++i) {
      sl_para_coord[i] = ScalarParametricCoordinate_{para_coord[i]};
    }

    int i{0};
    for (Index_ non_zero_basis_function{parameter_space.First()};
         non_zero_basis_function != parameter_space.Behind();
         ++non_zero_basis_function) {

      Index_ const& basis_function =
          (parameter_space.FindFirstNonZeroBasisFunction(sl_para_coord)
           + non_zero_basis_function.GetIndex());

      const auto evaluated =
          parameter_space.EvaluateBasisFunction(basis_function,
                                                non_zero_basis_function,
                                                sl_para_coord,
                                                unique_evaluations);

      basis[i] = evaluated;
      support[i] = basis_function.GetIndex1d().Get();
      ++i;
    }
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

  virtual void SplinepyInsertKnot(const int& p_dim, const double& knot) {
    splinepy::splines::helpers::ScalarTypeInsertKnot(*this, p_dim, knot);
  }

  virtual bool SplinepyRemoveKnot(const int& p_dim,
                                  const double& knot,
                                  const double& tolerance) {
    return splinepy::splines::helpers::ScalarTypeRemoveKnot(*this,
                                                            p_dim,
                                                            knot,
                                                            tolerance);
  }

  /// Bezier patch extraction
  virtual std::vector<std::shared_ptr<SplinepyBase>>
  SplinepyExtractBezierPatches() const {
    return splinepy::splines::helpers::ExtractBezierPatches<true>(*this);
  }

  constexpr const ParameterSpace_& GetParameterSpace() const {
    return *Base_::Base_::parameter_space_;
  }

  constexpr const VectorSpace_& GetVectorSpace() const {
    return *Base_::vector_space_;
  }

  constexpr Proximity_& GetProximity() { return *proximity_; }
  constexpr const Proximity_& GetProximity() const { return *proximity_; }

protected:
  std::unique_ptr<Proximity_> proximity_ = std::make_unique<Proximity_>(*this);
};

} /* namespace splinepy::splines */
