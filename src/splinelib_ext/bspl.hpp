#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

// SplineLib
#include <Sources/Splines/b_spline.hpp>

namespace py = pybind11;

template<int para_dim, int dim>
class BSplineExt : public splinelib::sources::splines::BSpline<para_dim, dim> {
public:
  using Base_ = splinelib::sources::splines::BSpline<para_dim, dim>;
  using Coordinate_ = typename Base_::Coordinate_;
  using Derivative_ = typename Base_::Derivative_;
  using Knot_ = typename Base_::Knot_;
  using ParameterSpace_ = typename Base_::ParameterSpace_;
  using ParametricCoordinate_ = typename Base_::ParametricCoordinate_;
  using VectorSpace_ = typename Base_::VectorSpace_;
  using OutputInformation_ = splinelib::Tuple<
      typename ParameterSpace_::OutputInformation_,
      typename VectorSpace_::OutputInformation_
  >;

  // Some private ones. Here we make it public
  using BezierInformation_ =
      typename Base_::ParameterSpace_::BezierInformation_;
  using BinomialRatios_ = typename Base_::ParameterSpace_::BinomialRatios_;
  using Index_ = typename Base_::Base_::Index_;
  using IndexLength_ = typename Index_::Length_;
  using IndexValue_ = typename Index_::Value_;
  using KnotRatios_ = typename Base_::ParameterSpace_::KnotRatios_;
  using Knots_ = typename Base_::Base_::Knots_;
  using BinomialRatio_ = typename BinomialRatios_::value_type;
  using KnotRatio_ = typename KnotRatios_::value_type;

  // Constructor
  using Base_::Base_;

  // update degrees since its size never changes
  void UpdateDegrees(int* p_degree_ptr) {
    ParameterSpace_ const &parameter_space = *Base_::Base_::parameter_space_;
    for (int i = 0; i < para_dim; i++) {
      p_degree_ptr[i] = parameter_space.degrees_[i].Get(); 
    }
  }

  // update current knot vectors to python
  // since list is mutable, update works
  void UpdateKnotVectors(py::list &p_knot_vectors) {

    // start clean
    p_knot_vectors.attr("clear")();

    ParameterSpace_ const &parameter_space = *Base_::Base_::parameter_space_;
    for (auto& knotvector : parameter_space.knot_vectors_) {
      auto const &kv = *knotvector; // in
      py::list p_kv; // out
      for (int i = 0; i < kv.GetSize(); i++) {
        auto const &knot = kv[splinelib::Index{i}];
        p_kv.append(knot.Get());
      }
      p_knot_vectors.append(p_kv);
    }
  }

  int GetNCps() {
    VectorSpace_ const &vector_space = *Base_::vector_space_;
    return vector_space.GetNumberOfCoordinates();
  }

  // control point changes in size, but you can figure it out with GetNCps()
  void UpdateControlPoints(double* cps_buf_ptr) {
    VectorSpace_ const &vector_space = *Base_::vector_space_;
    int numcps = vector_space.GetNumberOfCoordinates();

    // fill it up, phil!
    for (int i = 0; i < numcps; i++) {
      auto const &coord_named_phil = vector_space[splinelib::Index{i}];
      for (int j = 0; j < dim; j++) {
        cps_buf_ptr[i * dim + j] = coord_named_phil[j].Get();
      }
    }
  }

  // Computes (degree + 1) X ...
  // adapted from `SplineLib/Sources/Splines/b_spline.inc` 
  void BasisFunctionsAndIDs(ParametricCoordinate_ const &parametric_coordinate,
                            double* basis_function_values,
                            int* support_control_point_ids) const {

    ParameterSpace_ const &parameter_space = *Base_::Base_::parameter_space_;

    int i = 0;
    for (Index_ non_zero_basis_function{parameter_space.First()};
         non_zero_basis_function != parameter_space.Behind();
         ++non_zero_basis_function)
    {

      Index_ const &basis_function = (
          parameter_space.FindFirstNonZeroBasisFunction(parametric_coordinate)
          + non_zero_basis_function.GetIndex()
      );

      const auto evaluated = parameter_space.EvaluateBasisFunction(
          basis_function,
          parametric_coordinate
      );

      basis_function_values[i] = evaluated;
      support_control_point_ids[i] = basis_function.GetIndex1d().Get();
      i++;
    }
  }


};

