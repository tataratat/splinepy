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

    //c_bezier.UpdateDegrees(c_degrees);

    // update cps
    py::buffer_info cps_buf = p_control_points.request();
    double* cps_buf_ptr = static_cast<double *>(cps_buf.ptr);
    for (int i = 0; i < cps_buf.shape[0]; i++) {
      for (int j = 0; j < dim; j++) {
        c_bezier.control_points[i][j] = cps_buf_ptr[i * dim + j];
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
    for (int i = 0; i < cps_buf.shape[0]; i++) {
      for (int j = 0; j < dim; j++) {
        cps_buf_ptr[i * dim + j] = c_bezier.control_points[i][j];
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
      beziermanipulation::Point<static_cast<unsigned>(para_dim)> qpt; // query
      for (int j = 0; j < para_dim; j++) {
         qpt[j] = q_buf_ptr[i * para_dim + j];
      }
      const auto& eqpt = c_bezier.ForwardEvaluate(qpt); // evaluated query pt
      for (int k = 0; k < dim; k++) {
        r_buf_ptr[i * dim + k] = eqpt[k];
      }
    }

    results.resize({(int) q_buf.shape[0], dim});

    return results;
  }

  py::array_t<double> classique_evaluate(py::array_t<double> queries) {
    update_unless_skip();

    py::buffer_info q_buf = queries.request();
    double* q_buf_ptr = static_cast<double *>(q_buf.ptr);

    // prepare result arr
    py::array_t<double> results(q_buf.shape[0] * dim);
    py::buffer_info r_buf = results.request();
    double* r_buf_ptr = static_cast<double *>(r_buf.ptr);

    for (int i = 0; i < q_buf.shape[0]; i++) {
      beziermanipulation::Point<static_cast<unsigned>(para_dim)> qpt; // query
      for (int j = 0; j < para_dim; j++) {
         qpt[j] = q_buf_ptr[i * para_dim + j];
      }
      const auto& eqpt = c_bezier.Evaluate(qpt); // evaluated query pt
      for (int k = 0; k < dim; k++) {
        r_buf_ptr[i * dim + k] = eqpt[k];
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
        .def("classique_evaluate",
                 &PyBezier<para_dim, dim>::classique_evaluate,
                 py::arg("queries"))
        .def("elevate_degree",
                 &PyBezier<para_dim, dim>::elevate_degree,
                 py::arg("p_dim"))
        ;
}
