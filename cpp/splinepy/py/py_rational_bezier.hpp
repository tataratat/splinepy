#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <bezman/src/bezier_group.hpp>
#include <bezman/src/point.hpp>
#include <bezman/src/rational_bezier_spline.hpp>
#include <splinepy/py/py_bezier.hpp>
#include <type_traits>

namespace py = pybind11;

template<std::size_t para_dim, std::size_t dim>
using RationalBezier = bezman::RationalBezierSpline<
    static_cast<std::size_t>(para_dim),
    std::conditional_t<(dim > 1),
                       bezman::Point<static_cast<unsigned>(dim)>,
                       double>,
    double>;

template<std::size_t para_dim, std::size_t dim>
class PyRationalBezier {
private:
  // Alias to the internal RationalBezier type
  using RationalBezierSpline_ = RationalBezier<para_dim, dim>;

public:
  int para_dim_ = para_dim;
  int dim_ = dim;

  // python arrays
  py::array_t<int> p_degrees;
  py::array_t<double> p_control_points;
  py::array_t<double> p_weights;

  /// Updating flag
  bool skip_update = false;

  /// C++ object from bezman backend
  RationalBezierSpline_ c_rational_bezier;

  /// Identifier
  const std::string whatami =
      "RationalBezier, parametric dimension: " + std::to_string(para_dim)
      + ", physical dimension: " + std::to_string(dim);

  /// Empty constructor
  PyRationalBezier() = default;

  /// Initializing constructor
  PyRationalBezier(
      // degrees
      py::array_t<int> degrees,
      // control points
      py::array_t<double> control_points,
      // weights
      py::array_t<double> weights)
      : p_degrees{degrees},
        p_control_points{control_points},
        p_weights{weights} {
    update_c();
  }

  /// Copy constructor - (can be defaulted?)
  PyRationalBezier(RationalBezierSpline_ rhs) {
    // Init c_rational_bezier using move constructor
    c_rational_bezier = std::move(rhs);
    p_control_points.resize(
        {(int) c_rational_bezier.GetWeightedControlPoints().size(), (int) dim});
    p_weights.resize({(int) c_rational_bezier.GetWeights().size(), 1});
    p_degrees.resize({(int) para_dim});
    update_p();
  }

  /// Transformation Constructor
  PyRationalBezier(const PyBezier<para_dim, dim>& rhs)
      : c_rational_bezier{rhs.c_bezier} {
    p_control_points.resize(
        {(int) c_rational_bezier.GetWeightedControlPoints().size(), (int) dim});
    p_weights.resize({(int) c_rational_bezier.GetWeights().size(), 1});
    p_degrees.resize({(int) para_dim});
    update_p();
  }

  /// Update Backend-object
  void update_c() {
    // Update the degrees from python -> c++
    std::array<std::size_t, para_dim> c_degrees;
    int* ds_ptr = static_cast<int*>(p_degrees.request().ptr);
    std::copy_n(ds_ptr, para_dim, c_degrees.begin());

    // (re)init c++ object
    c_rational_bezier = std::move(RationalBezierSpline_{c_degrees});

    // update ctps and weights
    double* cps_ptr = static_cast<double*>(p_control_points.request().ptr);
    double* weights_ptr = static_cast<double*>(p_weights.request().ptr);
    // Check if size matches
    assert(p_control_points.request().shape[0]
           == c_rational_bezier.NumberOfControlPoints);
    for (std::size_t i = 0;
         i < static_cast<std::size_t>(p_control_points.request().shape[0]);
         i++) {
      // Update weight
      c_rational_bezier.GetWeights()[i] = weights_ptr[i];
      // Update CTPS
      if constexpr (dim > 1) {
        for (std::size_t j = 0; j < dim; j++) {
          c_rational_bezier.GetWeightedControlPoints()[i][j] =
              cps_ptr[i * dim + j] * weights_ptr[i];
        }
      } else {
        c_rational_bezier.GetWeightedControlPoints()[i] =
            cps_ptr[i] * weights_ptr[i];
      }
    }
  }

  /// Update Python side
  void update_p() {
    // degrees
    int* ds_ptr = static_cast<int*>(p_degrees.request().ptr);
    // control points
    double* cps_ptr = static_cast<double*>(p_control_points.request().ptr);
    // weights
    double* weights_ptr = static_cast<double*>(p_weights.request().ptr);

    // update degrees
    for (std::size_t i = 0; i < para_dim; i++) {
      ds_ptr[i] = c_rational_bezier.GetDegrees()[i];
    }

    // update control_points
    // Check if shape changed
    if (static_cast<long int>(
            c_rational_bezier.GetWeightedControlPoints().size())
        != p_control_points.request().shape[0]) {
      const std::size_t number_of_ctps =
          c_rational_bezier.GetWeightedControlPoints().size();
      // Update Control Point Vector
      p_control_points = py::array_t<double>(number_of_ctps * dim);
      p_control_points.resize({(int) number_of_ctps, (int) dim});
      // Update Control Point Vector
      p_weights = py::array_t<double>(number_of_ctps);
      p_weights.resize({(int) number_of_ctps, 1});
      // Update pointers
      cps_ptr = static_cast<double*>(p_control_points.request().ptr);
      weights_ptr = static_cast<double*>(p_weights.request().ptr);
    }

    // update control_points
    for (std::size_t i = 0;
         i < static_cast<std::size_t>(p_control_points.request().shape[0]);
         i++) {
      weights_ptr[i] = c_rational_bezier.GetWeights()[i];
      double inv_weight_i = static_cast<double>(1.) / weights_ptr[i];
      if constexpr (dim > 1) {
        for (std::size_t j = 0; j < dim; j++) {
          cps_ptr[i * dim + j] =
              c_rational_bezier.GetWeightedControlPoints()[i][j] * inv_weight_i;
        }
      } else {
        cps_ptr[i] =
            c_rational_bezier.GetWeightedControlPoints()[i] * inv_weight_i;
      }
    }
  }

  /// Conditional update
  void update_unless_skip() {
    if (!skip_update) {
      update_c();
    }
  }

  // Evaluate
  py::array_t<double> evaluate(py::array_t<double> queries) {
    update_unless_skip();

    // Access query points
    double* q_ptr = static_cast<double*>(queries.request().ptr);

    // initialize results
    py::array_t<double> results(queries.request().shape[0] * dim);
    double* r_ptr = static_cast<double*>(results.request().ptr);

    // Loop over query points for evaluation
    for (std::size_t i = 0;
         i < static_cast<std::size_t>(queries.request().shape[0]);
         i++) {
      bezman::Point<static_cast<unsigned>(para_dim)> query_ptr;
      for (std::size_t j = 0; j < para_dim; j++) {
        query_ptr[j] = q_ptr[i * para_dim + j];
      }
      const auto& eqpt = c_rational_bezier.ForwardEvaluate(query_ptr);
      if constexpr (dim > 1) {
        for (std::size_t j = 0; j < dim; j++) {
          r_ptr[i * dim + j] = eqpt[j];
        }
      } else {
        r_ptr[i * dim] = eqpt;
      }
    }

    results.resize({(int) queries.request().shape[0], (int) dim});

    return results;
  }

  // Elevate the order along one specific parametric dimension
  void elevate_degree(int p_dim) {
    update_unless_skip();

    c_rational_bezier.OrderElevateAlongParametricDimension(
        static_cast<std::size_t>(p_dim));

    update_p();
  }

  // Multiplication routines
  PyRationalBezier<para_dim, 1>
  multiply_with_spline(const PyRationalBezier& a) {
    PyRationalBezier<para_dim, 1> result{(*this).c_rational_bezier
                                         * a.c_rational_bezier};
    result.update_p();
    return result;
  }

  // Addition routines
  PyRationalBezier add_spline(const PyRationalBezier& a) {
    PyRationalBezier result{(*this).c_rational_bezier + a.c_rational_bezier};
    result.update_p();
    return result;
  }

  PyRationalBezier
  multiply_with_scalar_spline(const PyRationalBezier<para_dim, 1>& a) {
    PyRationalBezier result{(*this).c_rational_bezier * a.c_rational_bezier};
    result.update_p();
    return result;
  }

  // Addition Routines

  // Composition Routine
  template<int par_dim_inner_function>
  PyRationalBezier<par_dim_inner_function, dim>
  ComposeRR(const PyRationalBezier<par_dim_inner_function, para_dim>&
                inner_function) {
    // Use Composition routine
    PyRationalBezier<par_dim_inner_function, dim> result{
        (*this).c_rational_bezier.Compose(inner_function.c_rational_bezier)};
    result.update_p();
    return result;
  }

  // Composition Routine
  template<int par_dim_inner_function>
  PyRationalBezier<par_dim_inner_function, dim>
  ComposeRP(const PyBezier<par_dim_inner_function, para_dim>& inner_function) {
    // Use Composition routine
    PyRationalBezier<par_dim_inner_function, dim> result{
        (*this).c_rational_bezier.Compose(inner_function.c_bezier)};
    result.update_p();
    return result;
  }
};

// Start defining the python interface
template<int para_dim, int dim>
void add_rational_bezier_pyclass(py::module& m, const char* class_name) {
  py::class_<PyRationalBezier<para_dim, dim>> klasse(m, class_name);

  klasse
      // Default initializer
      .def(py::init<>())
      // initializer with values
      .def(py::init<py::array_t<int>,
                    py::array_t<double>,
                    py::array_t<double>>(),
           py::arg("degrees"),
           py::arg("control_points"),
           py::arg("weights"))
      // Transformation Copy constructor
      .def(py::init<PyBezier<para_dim, dim>>())
      // Identifier
      .def_readonly("whatami", &PyRationalBezier<para_dim, dim>::whatami)
      // physical dimension
      .def_readonly("dim", &PyRationalBezier<para_dim, dim>::dim_)
      // Parametric dimension
      .def_readonly("para_dim", &PyRationalBezier<para_dim, dim>::para_dim_)
      // Skip update
      .def_readwrite("skip_update",
                     &PyRationalBezier<para_dim, dim>::skip_update)
      // Degrees
      .def_readwrite("degrees", &PyRationalBezier<para_dim, dim>::p_degrees)
      // Control Point Vector
      .def_readwrite("control_points",
                     &PyRationalBezier<para_dim, dim>::p_control_points)
      // RWeights Vector
      .def_readwrite("weights", &PyRationalBezier<para_dim, dim>::p_weights)
      // Update Backend data
      .def("update_c", &PyRationalBezier<para_dim, dim>::update_c)
      // Update Python side
      .def("update_p", &PyRationalBezier<para_dim, dim>::update_p)
      // Evaluate
      .def("evaluate",
           &PyRationalBezier<para_dim, dim>::evaluate,
           py::arg("queries"))
      // Degree elevation
      .def("elevate_degree",
           &PyRationalBezier<para_dim, dim>::elevate_degree,
           py::arg("p_dim"))
      // Multiplication
      .def("multiply_with_spline",
           &PyRationalBezier<para_dim, dim>::multiply_with_spline,
           py::arg("factor"))
      // Multiplication with a scalar spline
      .def("multiply_with_scalar_spline",
           &PyRationalBezier<para_dim, dim>::multiply_with_scalar_spline,
           py::arg("factor"))
      // Addition
      .def("add_spline",
           &PyRationalBezier<para_dim, dim>::add_spline,
           py::arg("summand"))
      // Composition with one dimensional spline
      .def("compose_line_rr",
           &PyRationalBezier<para_dim, dim>::template ComposeRR<1>,
           py::arg("inner_function"))
      // Composition with a surface
      .def("compose_surface_rr",
           &PyRationalBezier<para_dim, dim>::template ComposeRR<2>,
           py::arg("inner_function"))
      // Composition with a volume
      .def("compose_volume_rr",
           &PyRationalBezier<para_dim, dim>::template ComposeRR<3>,
           py::arg("inner_function"))
      // Composition with one dimensional spline
      .def("compose_line_rp",
           &PyRationalBezier<para_dim, dim>::template ComposeRP<1>,
           py::arg("inner_function"))
      // Composition with a surface
      .def("compose_surface_rp",
           &PyRationalBezier<para_dim, dim>::template ComposeRP<2>,
           py::arg("inner_function"))
      // Composition with a volume
      .def("compose_volume_rp",
           &PyRationalBezier<para_dim, dim>::template ComposeRP<3>,
           py::arg("inner_function"))
      // Picke
      .def(py::pickle(
          [](const PyRationalBezier<para_dim, dim>& RationalBezier) {
            return py::make_tuple(RationalBezier.p_degrees,
                                  RationalBezier.p_control_points,
                                  RationalBezier.p_weights);
          },
          [](py::tuple t) {
            if (t.size() != 3) {
              throw std::runtime_error("Invalid PyRationalBezier state!");
            }

            PyRationalBezier<para_dim, dim> pyb(
                t[0].cast<py::array_t<int>>(),
                t[1].cast<py::array_t<double>>(),
                t[2].cast<py::array_t<double>>());

            return pyb;
          }));
}
