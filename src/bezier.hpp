#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <bezierManipulation/src/point.hpp>
#include <bezierManipulation/src/bezier_spline.hpp>

namespace py = pybind11;

template<int para_dim, int dim>
using Bezier = beziermanipulation::BezierSpline<
    static_cast<std::size_t>(para_dim),
    beziermanipulation::Point<static_cast<unsigned>(dim)>,
    double
>;


template<int para_dim, int dim>
class PyBezier {
public:

  //using ScalarType_ = Bezier<para_dim, dim>;

  // python arrays
  py::array_t<int> p_degrees;
  py::array_t<double> p_control_points;

  bool skip_update = false;

  Bezier<para_dim, dim> c_spline;

  const std::string whatami =
      "Bezier, parametric dimension: "
      + std::to_string(para_dim)
      + ", physical dimension: "
      + std::to_string(dim);

  PyBezier();
  PyBezier(py::array_t<int> degrees,
           py::array_t<double> control_points) :
                p_degrees(degrees),
                p_control_points(control_points) {
    update_c();
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
    c_spline = std::move(Bezier<para_dim, dim>{c_degrees});

    // update cps
    py::buffer_info cps_buf = p_control_points.request();
    double* cps_buf_ptr = static_cast<double *>(cps_buf.ptr);
    for (int i = 0; i < cps_buf.shape[0]; i++) {
      for (int j = 0; j < dim; j++) {
        c_spline.control_points[i][j] = cps_buf_ptr[i * dim + j];
      }
    }

  }

  void update_p() {
    py::buffer_info ds_buf = p_degrees.request();
    int* ds_buf_ptr = static_cast<int *>(ds_buf.ptr);
   
  }

  py::array_t<double> evaluate(py::array_t<double> queries) {
    py::buffer_info q_buf = queries.request();
    double* q_buf_ptr = static_cast<double *>(q_buf.ptr);

    

    for (int i = 0; i < q_buf.shape[0]; i++) {
      //for (int j = 0; 
    }
  }

};

template<int para_dim, int dim>
void add_bezier_pyclass(py::module &m, const char *class_name) {
  py::class_<PyBezier<para_dim, dim>> klasse(m, class_name);

  klasse.def(py::init<>())
        //.def(py::init<py::array_t<int>, py::array_t<double>>(),
        //         py::arg("degrees"),
        //        py::arg("control_points"))
        ;
}
