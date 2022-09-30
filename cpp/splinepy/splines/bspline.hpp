#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

// SplineLib
#include <Sources/Splines/b_spline.hpp>

#include <splinepy/proximity/proximity.hpp>
#include <splinepy/splines/helpers/properties.hpp>
#include <splinepy/splines/helpers/scalar_type_wrapper.hpp>
#include <splinepy/splines/splinepy_base.hpp>

namespace splinepy::splines {

namespace py = pybind11;

template<int para_dim, int dim>
class BSpline : public splinepy::splines::SplinepyBase,
                public splinelib::sources::splines::BSpline<para_dim, dim> {
public:
  static constexpr int kParaDim = para_dim;
  static constexpr int kDim = dim;

  // TODO remve afterwards
  constexpr static int para_dim_ = para_dim;
  constexpr static int dim_ = dim;

  // splinepy
  using SplinepyBase_ = typename splinepy::splines::SplinepyBase;

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
  using Dimension_ = typename splinelib::Dimension;
  using Tolerance_ = typename splinelib::sources::splines::Tolerance;
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

  // raw ptr based inithelper.
  // degrees should have same size as parametric dimension
  // having knot_vectors vector of vector, we can keep track of their length,
  // as well as the legnth of control_points/weights.
  Base_ RawPtrInitHelper(const double* degrees,
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
  BSpline(const double* degrees,
          const std::vector<std::vector<double>>& knot_vectors,
          const double* control_points)
      : Base_(RawPtrInitHelper(degrees, knot_vectors, control_points)) {}
  // inherited ctor
  using Base_::Base_;

  // required implementations
  virtual constexpr int SplinepyParaDim() const { return kParaDim; }

  virtual constexpr int SplinepyDim() const { return kDim; }

  virtual std::string SplinepyWhatAmI() const {
    return "BSpline, parametric dimension: " + std::to_string(SplinepyParaDim())
           + ", physical dimension: " + std::to_string(SplinepyDim());
  }

  virtual int SplinepyNumberOfControlPoints() const {
    return GetVectorSpace().GetNumberOfCoordinates();
  }

  virtual void
  SplinepyCurrentProperties(double* degrees,
                            std::vector<std::vector<double>>* knot_vectors,
                            double* control_points,
                            double* weights /* untouched */) const {

    const auto& parameter_space = GetParameterSpace();
    const auto& vector_space = GetVectorSpace();

    // degrees
    for (std::size_t i{}; i < kParaDim; ++i) {
      degrees[i] = static_cast<double>(parameter_space.GetDegrees()[i]);
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
      for (int j{}; j < kvsize; ++j) {
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

  const ParameterSpace_& GetParameterSpace() const {
    return *Base_::Base_::parameter_space_;
  }

  const VectorSpace_& GetVectorSpace() const { return *Base_::vector_space_; }

  // update degrees since its size never changes
  void UpdateDegrees(int* p_degree_ptr) {
    ParameterSpace_ const& parameter_space = *Base_::Base_::parameter_space_;
    for (int i = 0; i < para_dim; i++) {
      p_degree_ptr[i] = parameter_space.GetDegrees()[i].Get();
    }
  }

  // update current knot vectors to python
  // since list is mutable, update works
  void UpdateKnotVectors(py::list& p_knot_vectors) {

    // start clean
    p_knot_vectors.attr("clear")();

    ParameterSpace_ const& parameter_space = *Base_::Base_::parameter_space_;
    for (auto& knotvector : parameter_space.GetKnotVectors()) {
      auto const& kv = *knotvector; // in
      py::list p_kv;                // out
      for (int i = 0; i < kv.GetSize(); i++) {
        auto const& knot = kv[splinelib::Index{i}];
        p_kv.append(knot.Get());
      }
      p_knot_vectors.append(p_kv);
    }
  }

  int GetNCps() {
    VectorSpace_ const& vector_space = *Base_::vector_space_;
    return vector_space.GetNumberOfCoordinates();
  }

  // control point changes in size, but you can figure it out with GetNCps()
  void UpdateControlPoints(double* cps_buf_ptr) {
    VectorSpace_ const& vector_space = *Base_::vector_space_;
    int numcps = vector_space.GetNumberOfCoordinates();

    // fill it up, phil!
    for (int i = 0; i < numcps; i++) {
      auto const& coord_named_phil = vector_space[splinelib::Index{i}];
      for (int j = 0; j < dim; j++) {
        cps_buf_ptr[i * dim + j] = coord_named_phil[j].Get();
      }
    }
  }

  // control point changes in size, but you can figure it out with GetNCps()
  void FillControlPoints(double* cps_buf_ptr) const {
    VectorSpace_ const& vector_space = *Base_::vector_space_;
    int numcps = vector_space.GetNumberOfCoordinates();

    // fill it up, phil!
    for (int i = 0; i < numcps; i++) {
      auto const& coord_named_phil = vector_space[splinelib::Index{i}];
      for (int j = 0; j < dim; j++) {
        cps_buf_ptr[i * dim + j] = coord_named_phil[j].Get();
      }
    }
  }

  // Computes (degree + 1) X ...
  // adapted from `SplineLib/Sources/Splines/b_spline.inc`
  void BasisFunctionsAndIDs(ParametricCoordinate_ const& parametric_coordinate,
                            double* basis_function_values,
                            int* support_control_point_ids) const {

    ParameterSpace_ const& parameter_space = *Base_::Base_::parameter_space_;

    int i = 0;
    for (Index_ non_zero_basis_function{parameter_space.First()};
         non_zero_basis_function != parameter_space.Behind();
         ++non_zero_basis_function) {

      Index_ const& basis_function =
          (parameter_space.FindFirstNonZeroBasisFunction(parametric_coordinate)
           + non_zero_basis_function.GetIndex());

      const auto evaluated =
          parameter_space.EvaluateBasisFunction(basis_function,
                                                parametric_coordinate);

      basis_function_values[i] = evaluated;
      support_control_point_ids[i] = basis_function.GetIndex1d().Get();
      i++;
    }
  }

  Proximity_& GetProximity() {
    if (!proximity_initialized_) {
      proximity_ = std::make_unique<Proximity_>(*this);
      proximity_initialized_ = true;
    }
    return *proximity_;
  }

protected:
  std::unique_ptr<Proximity_> proximity_;
  bool proximity_initialized_ = false;
};

/// dynamic creation of templated BSpline
std::shared_ptr<SplinepyBase> SplinepyBase::SplinepyCreateBSpline(
    const int para_dim,
    const int dim,
    const double* degrees,
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
    }
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
    }
  }
}

} /* namespace splinepy::splines */
