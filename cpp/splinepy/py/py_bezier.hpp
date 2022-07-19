#pragma once

#include<type_traits>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <bezman/src/point.hpp>
#include <bezman/src/bezier_spline.hpp>

namespace py = pybind11;

template<int para_dim, int dim>
using Bezier = bezman::BezierSpline<
    static_cast<std::size_t>(para_dim),
    std::conditional_t<
      (dim > 1),
      bezman::Point<static_cast<unsigned>(dim)>,
      double>,
    double
>;


template<int para_dim, int dim>
class PyBezier {
private:
  // Alias to the internal Bezier type
  using BezierSpline_ = Bezier<para_dim, dim>;

public:

  int para_dim_ = para_dim;
  int dim_ = dim;

  // python arrays
  py::array_t<int> p_degrees;
  py::array_t<double> p_control_points;

  bool skip_update = false;

  Bezier<para_dim, dim> c_bezier;

  const std::string whatami =
      "Bezier, parametric dimension: "
      + std::to_string(para_dim)
      + ", physical dimension: "
      + std::to_string(dim);

  PyBezier() = default;
  PyBezier(py::array_t<int> degrees,
           py::array_t<double> control_points) :
                p_degrees(degrees),
                p_control_points(control_points) {
    update_c();
  }

  PyBezier(BezierSpline_ rhs){
    // Init c_bezier using move constructor
    c_bezier = std::move(rhs);
    p_control_points.resize({(int) c_bezier.control_points.size(), dim});
    p_degrees.resize({(int) para_dim});
    update_p();
  }

  void update_c() {

    std::array<std::size_t, para_dim> c_degrees;

    // update degrees
    py::buffer_info ds_buf = p_degrees.request();
    int* ds_buf_ptr = static_cast<int *>(ds_buf.ptr);

    for (int i = 0; i < para_dim; i++) {
      c_degrees[i] = ds_buf_ptr[i];
    }
    
    // (re)init
    c_bezier = std::move(Bezier<para_dim, dim>{c_degrees});

    // update cps
    py::buffer_info cps_buf = p_control_points.request();
    double* cps_buf_ptr = static_cast<double *>(cps_buf.ptr);
    for (int i = 0; i < cps_buf.shape[0]; i++) {
      if constexpr (dim > 1) {
        for (int j = 0; j < dim; j++) {
          c_bezier.control_points[i][j] = cps_buf_ptr[i * dim + j];
        }
      } else {
          c_bezier.control_points[i] = cps_buf_ptr[i];
      }
    }

  }

  void update_p() {
    // degrees
    py::buffer_info ds_buf = p_degrees.request();
    int* ds_buf_ptr = static_cast<int *>(ds_buf.ptr);
    // control points
    py::buffer_info cps_buf = p_control_points.request();
    double* cps_buf_ptr = static_cast<double *>(cps_buf.ptr);

    // update degrees
    auto const& c_ds = c_bezier.GetDegrees();
    for (int i = 0; i < para_dim; i++) {
      ds_buf_ptr[i] = c_ds[i];
    }

    // update control_points
    // Check if shape changed
    if (c_bezier.control_points.size() != cps_buf.shape[0]) {
      p_control_points = py::array_t<double>(c_bezier.control_points.size() * dim);
      cps_buf = p_control_points.request();
      cps_buf_ptr = static_cast<double *>(cps_buf.ptr);
      p_control_points.resize({(int) c_bezier.control_points.size(), dim});
    }

    // Update point coordinates
    for (int i = 0; i < cps_buf.shape[0]; i++) {
      if constexpr (dim > 1){
        for (int j = 0; j < dim; j++) {
          cps_buf_ptr[i * dim + j] = c_bezier.control_points[i][j];
        }
      } else {
          cps_buf_ptr[i] = c_bezier.control_points[i];
      }
    }
  }

  void update_unless_skip() {
    if(!skip_update) {
      update_c();
    }
  }

  py::array_t<double> evaluate(py::array_t<double> queries) {
    update_unless_skip();

    py::buffer_info q_buf = queries.request();
    double* q_buf_ptr = static_cast<double *>(q_buf.ptr);

    // prepare result arr
    py::array_t<double> results(q_buf.shape[0] * dim);
    py::buffer_info r_buf = results.request();
    double* r_buf_ptr = static_cast<double *>(r_buf.ptr);

    for (int i = 0; i < q_buf.shape[0]; i++) {
      bezman::Point<static_cast<unsigned>(para_dim)> qpt; // query
      for (int j = 0; j < para_dim; j++) {
         qpt[j] = q_buf_ptr[i * para_dim + j];
      }
      const auto& eqpt = c_bezier.ForwardEvaluate(qpt); // evaluated query pt
      if constexpr (dim > 1){
        for (int j = 0; j < dim; j++) {
          r_buf_ptr[i * dim + j] = eqpt[j];
        }
      } else {
          r_buf_ptr[i * dim] = eqpt;
      }
    }

    results.resize({(int) q_buf.shape[0], dim});

    return results;
  }

  py::array_t<double> pseudorecursive_evaluate(py::array_t<double> queries) {
    update_unless_skip();

    py::buffer_info q_buf = queries.request();
    double* q_buf_ptr = static_cast<double *>(q_buf.ptr);

    // prepare result arr
    py::array_t<double> results(q_buf.shape[0] * dim);
    py::buffer_info r_buf = results.request();
    double* r_buf_ptr = static_cast<double *>(r_buf.ptr);

    for (int i = 0; i < q_buf.shape[0]; i++) {
      bezman::Point<static_cast<unsigned>(para_dim)> qpt; // query
      for (int j = 0; j < para_dim; j++) {
         qpt[j] = q_buf_ptr[i * para_dim + j];
      }
      const auto& eqpt = c_bezier.Evaluate(qpt); // evaluated query pt
      if constexpr (dim > 1){
        for (int j = 0; j < dim; j++) {
          r_buf_ptr[i * dim + j] = eqpt[j];
        }
      } else {
          r_buf_ptr[i * dim] = eqpt;
      }
    }

    results.resize({(int) q_buf.shape[0], dim});

    return results;
  }

  void elevate_degree(int p_dim) {
    update_unless_skip();

    c_bezier.OrderElevateAlongParametricDimension(
        static_cast<std::size_t>(p_dim)
    );

    update_p();
  }

  // Multiplication routines
  PyBezier<para_dim,1> multiply_with_spline (const PyBezier& a){
    PyBezier<para_dim,1> result{(*this).c_bezier * a.c_bezier};
    result.update_p();
    return result;
  }
  
  // Addition routines
  PyBezier add_spline (const PyBezier& a){
    PyBezier result{(*this).c_bezier + a.c_bezier};
    result.update_p();
    return result;
  }
  
  PyBezier multiply_with_scalar_spline (const PyBezier<para_dim,1>& a){
    PyBezier result{(*this).c_bezier * a.c_bezier};
    result.update_p();
    return result;
  }

  // Addition Routines

  // Composition Routine
  template<int par_dim_inner_function>
  PyBezier<par_dim_inner_function, dim> Compose(
              const PyBezier<par_dim_inner_function, para_dim>& inner_function){
      // Use Composition routine
      PyBezier<par_dim_inner_function, dim> result{(*this).c_bezier.Compose(
                        inner_function.c_bezier)
                        };
    result.update_p();
    return result;
  }


};

template<int para_dim, int dim>
void add_bezier_pyclass(py::module &m, const char *class_name) {
  py::class_<PyBezier<para_dim, dim>> klasse(m, class_name);

  klasse.def(py::init<>())
        .def(py::init<py::array_t<int>, py::array_t<double>>(),
                 py::arg("degrees"),
                 py::arg("control_points"))
        .def_readonly("whatami",
                           &PyBezier<para_dim, dim>::whatami)
        .def_readonly("dim",
                           &PyBezier<para_dim, dim>::dim_)
        .def_readonly("para_dim",
                           &PyBezier<para_dim, dim>::para_dim_)
        .def_readwrite("skip_update",
                           &PyBezier<para_dim, dim>::skip_update)
        .def_readwrite("degrees",
                           &PyBezier<para_dim, dim>::p_degrees)
        .def_readwrite("control_points",
                           &PyBezier<para_dim, dim>::p_control_points)
        .def("update_c",
                 &PyBezier<para_dim, dim>::update_c)
        .def("update_p",
                 &PyBezier<para_dim, dim>::update_p)
        .def("evaluate",
                 &PyBezier<para_dim, dim>::evaluate,
                 py::arg("queries"))
        .def("recursive_evaluate",
                 &PyBezier<para_dim, dim>::pseudorecursive_evaluate,
                 py::arg("queries"))
        .def("elevate_degree",
                 &PyBezier<para_dim, dim>::elevate_degree,
                 py::arg("p_dim"))
        .def("multiply_with_spline",
                 &PyBezier<para_dim, dim>::multiply_with_spline,
                 py::arg("factor"))
        .def("multiply_with_scalar_spline",
                 &PyBezier<para_dim, dim>::multiply_with_scalar_spline,
                 py::arg("factor"))
        .def("add_spline",
                 &PyBezier<para_dim, dim>::add_spline,
                 py::arg("summand"))
        .def("compose_line",
                 &PyBezier<para_dim, dim>::template Compose<1>,
                 py::arg("inner_function"))
        .def("compose_surface",
                 &PyBezier<para_dim, dim>::template Compose<2>,
                 py::arg("inner_function"))
        .def("compose_volume",
                 &PyBezier<para_dim, dim>::template Compose<3>,
                 py::arg("inner_function"))
        .def(py::pickle(
                 [] (const PyBezier<para_dim, dim> &bezier) {
                   return py::make_tuple(
                       bezier.p_degrees,
                       bezier.p_control_points
                   );
                 },
                 [] (py::tuple t) {
                   if (t.size() != 2) {
                     throw std::runtime_error("Invalid PyBezier state!");
                   }

                   PyBezier<para_dim, dim> pyb(
                     t[0].cast<py::array_t<int>>(),
                     t[1].cast<py::array_t<double>>()
                   );

                   return pyb;
                 }
             ));
        
}
