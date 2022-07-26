#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

// SplineLib
#include <Sources/Splines/nurbs.hpp>

#include <splinepy/proximity/proximity.hpp>

namespace splinepy::splines {

namespace py = pybind11;

template<int para_dim, int dim>
class Nurbs : public splinelib::sources::splines::Nurbs<para_dim, dim> {
public:
  static constexpr int kParaDim = para_dim;
  static constexpr int kDim = dim;

  using Base_ = splinelib::sources::splines::Nurbs<para_dim, dim>;
  using Coordinate_ = typename Base_::Coordinate_;
  using ScalarCoordinate_ = typename Coordinate_::value_type;
  using Derivative_ = typename Base_::Derivative_;
  using Knot_ = typename Base_::Knot_;
  using ParameterSpace_ = typename Base_::ParameterSpace_;
  using ParametricCoordinate_ = typename Base_::ParametricCoordinate_;
  using ScalarParametricCoordinate_ =
      typename ParametricCoordinate_::value_type;
  using WeightedVectorSpace_ = typename Base_::WeightedVectorSpace_;
  using OutputInformation_ = splinelib::Tuple<
      typename ParameterSpace_::OutputInformation_,
      typename WeightedVectorSpace_::OutputInformation_
  >;
  using VectorSpace_ = typename WeightedVectorSpace_::Base_; // <dim + 1>

  // Some private ones. Here we make it public
  using Index_ = typename Base_::Base_::Index_;
  using IndexLength_ = typename Index_::Length_;
  using IndexValue_ = typename Index_::Value_;
  using KnotRatios_ = typename Base_::ParameterSpace_::KnotRatios_;
  using Knots_ = typename Base_::Base_::Knots_;
  using HomogeneousBSpline_ = typename Base_::HomogeneousBSpline_;

  using Proximity_ = splinepy::proximity::Proximity<Nurbs<para_dim, dim>>;

  // Constructor
  using Base_::Base_;

  constexpr static int para_dim_ = para_dim;
  constexpr static int dim_ = dim;

  const ParameterSpace_& GetParameterSpace() const {
    return *Base_::Base_::parameter_space_;
  }

  const WeightedVectorSpace_& GetWeightedVectorSpace() const {
    return *Base_::weighted_vector_space_;
  }

   /* 
    * Update degrees since its size never changes
   *
   * Parameters
   * -----------
   * ds_buf_ptr: out
   */
  void UpdateDegrees(int* ds_buf_ptr) {

    ParameterSpace_ const &parameter_space = *Base_::Base_::parameter_space_;
    for (int i = 0; i < para_dim; i++) {
      ds_buf_ptr[i] = parameter_space.GetDegrees()[i].Get();
    }
  }

  /*
   * Update current knot vectors to python 
   * since list is mutable, update works
   *
   * Parameters
   * -----------
   * p_knot_vectors: inout
   *  will be `clear`ed before writing.
   */
  void UpdateKnotVectors(py::list &p_knot_vectors) {

    // start clean
    p_knot_vectors.attr("clear")();

    ParameterSpace_ const &parameter_space = *Base_::Base_::parameter_space_;
    for (auto& knotvector : parameter_space.GetKnotVectors()) {
      auto const &kv = *knotvector; // in
      py::list p_kv; // out
      for (int i = 0; i < kv.GetSize(); i++) {
        auto const &knot = kv[splinelib::Index{i}];
        p_kv.append(knot.Get());
      }
      p_knot_vectors.append(p_kv);
    }

  }

  /* Returns n_cpts
   * used to allocate numpy array
   *
   * Parameters
   * -----------
   * None
   *
   * Returns
   * -------
   * N_CPS
   */ 
  int GetNCps() {
    WeightedVectorSpace_ const &vector_space = *Base_::weighted_vector_space_;
    return vector_space.GetNumberOfCoordinates();
  }

  /* Update cps and weights at the same time since they belong together in
   * weightedvectorspace
   *
   * Parameters
   * -----------
   * cps_buf_ptr: out
   *   size -> (n_cps * dim)
   * ws_buf_ptr: out
   *   size -> (n_ws) (=n_cps) 
   */
  void UpdateControlPointsAndWeights(double* cps_buf_ptr,
                                     double* ws_buf_ptr) {
    WeightedVectorSpace_ const &vector_space = *Base_::weighted_vector_space_;
    int numcps = vector_space.GetNumberOfCoordinates();


    // fill it up, phil!
    for (int i = 0; i < numcps; i++) {
      auto const &coord_named_phil = vector_space[splinelib::Index{i}];
      // phil needs to first project before filling.
      auto const &projected_phil = WeightedVectorSpace_::Project(
          coord_named_phil
      );

      for (int j = 0; j < dim; j++) {
        cps_buf_ptr[i * dim + j] = projected_phil[j].Get();
      }
      ws_buf_ptr[i] = coord_named_phil[dim].Get();
    }

  }

  /* Given parametric coordinate, returns basis functions and support control
   * point ids.
   * pointer should have size of (degree + 1)^(para_dim)
   *
   * Parameters
   * -----------
   * parametric_coordinate: in
   * basis_function_values: out
   * support_control_point_ids: out
   */
  void BasisFunctionsAndIDs(ParametricCoordinate_ const &parametric_coordinate,
                            double* basis_function_values,
                            int* support_control_point_ids) const {

    ParameterSpace_ const &parameter_space = *Base_::Base_::parameter_space_;

    int i = 0;
    double W = 0.0;
    for (Index_ non_zero_basis_function{parameter_space.First()};
         non_zero_basis_function != parameter_space.Behind();
         ++non_zero_basis_function)
    {

      Index_ const &basis_function = (
          parameter_space.FindFirstNonZeroBasisFunction(parametric_coordinate) +
          non_zero_basis_function.GetIndex()
      );

      // general basis fn
      const auto evaluated = parameter_space.EvaluateBasisFunction(
          basis_function,
          parametric_coordinate
      );

      // get `w` and add to `W`
      const auto& support_id = basis_function.GetIndex1d(); // not int yet

      VectorSpace_ const &vector_space = *Base_::weighted_vector_space_;
      const auto& w = vector_space[support_id][dim].Get(); // dim: last elem

      const double N_times_w = evaluated * w;
 
      W += N_times_w;
      basis_function_values[i] = N_times_w; // not yet final
      support_control_point_ids[i] = support_id.Get();
      i++;
    }

    // Loop and divide entries by W
    int end = i;
    double W_inv = 1 / W;
    for (i = 0; i < end; i++) {
      basis_function_values[i] *= W_inv;
    }
  }

  /* Fn overload with pure array in,outs.
   *
   * Parameters
   * -----------
   * parametric_coordinate: in
   * basis_function_values: out
   * support_control_point_ids: out
   */
  void BasisFunctionsAndIDs(const double* parametric_coordinate,
                            double* basis_function_values,
                            int* support_control_point_ids) const {
    // prepare parametric coordinate
    ParametricCoordinate_ pc{};
    for (int i{0}; i < para_dim; i++) {
      pc[i] = ScalarParametricCoordinate_{parametric_coordinate[i]};
    }

    BasisFunctionsAndIDs(parametric_coordinate,
                         basis_function_values,
                         support_control_point_ids);
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

}; /* class Nurbs */

} /* namespace splinepy::splines */
