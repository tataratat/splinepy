#include <memory>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <tuple>
#include <string>

// pybind11
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

// SplineLib
#include <Sources/Splines/b_spline.hpp>
#include <Sources/Utilities/named_type.hpp>
#include <Sources/Utilities/string_operations.hpp>
#include <Sources/InputOutput/iges.hpp>
#include <Sources/InputOutput/xml.hpp>
#include <Sources/InputOutput/vtk.hpp>
#include <Sources/InputOutput/irit.hpp>

// Local
// Fitting
#include <fitting/fitting.hpp>

namespace py = pybind11;

using namespace splinelib::sources;

template<int para_dim, int dim>
struct PyBSpline {

    using BSpline = splines::BSpline<para_dim, dim>;

    // For writing cpp splines
    using ParameterSpace = typename BSpline::ParameterSpace_;
    using VectorSpace = typename BSpline::VectorSpace_;
    using Coordinates = typename VectorSpace::Coordinates_;
    using Degrees = typename ParameterSpace::Degrees_;
    using KnotVectors = typename ParameterSpace::KnotVectors_;
    using Coordinate = typename Coordinates::value_type;
    using Degree = typename Degrees::value_type;
    using KnotVector = typename KnotVectors::value_type::element_type;
    using Knots = typename KnotVector::Knots_;
    using ScalarCoordinate = typename Coordinate::value_type;
    using Knot = typename Knots::value_type;
    using KnotB = typename BSpline::Knot_;
    using ParametricCoordinate = typename BSpline::ParametricCoordinate_;
    using ScalarParametricCoordinate = typename ParametricCoordinate::value_type;
    using Derivative = typename BSpline::Derivative_;
    using NumberOfParametricCoordinates =  typename  ParameterSpace::NumberOfParametricCoordinates_;

    // For reading cpp splines
    using OutputInformation = typename BSpline::OutputInformation_;
    using OutputParameterSpace = typename std::tuple_element_t<0, OutputInformation>;
    using OutputVectorSpace = typename std::tuple_element_t<1, OutputInformation>;
    using OutputKnotVectors = typename std::tuple_element_t<0, OutputParameterSpace>;
    using OutputCoordinates = typename std::tuple_element_t<0, OutputVectorSpace>;
    using OutputDegrees = typename std::tuple_element_t<1, OutputParameterSpace>;

    // Counters
    int i = 0;
    int j = 0;

    // Hr. BSpline himself
    BSpline c_bspline;
    std::shared_ptr<ParameterSpace> c_parameter_space;
    std::shared_ptr<VectorSpace> c_vector_space;

    // Hr. BSpline's answer to the question "What am I?"
    const std::string whatami =
        "BSpline, parametric dimension: "
        + std::to_string(para_dim)
        + ", physical dimension: "
        + std::to_string(dim);

    // Hr. BSpline's python Family
    py::list p_knot_vectors;
    py::array_t<int> p_degrees;
    py::array_t<double> p_control_points;

    // Hr. BSpline's switch to skip_update_c
    bool skip_update = false;

    // Constructor.
    PyBSpline() {
        py::list empty_list;
        p_knot_vectors.attr("clear")();
        p_knot_vectors.append(empty_list);
    } 
    PyBSpline(py::array_t<int> degrees,
              py::list knot_vectors,
              py::array_t<double> control_points) :
                                                        p_degrees(degrees),
                                              p_knot_vectors(knot_vectors),
                                          p_control_points(control_points) {
        update_c();
    } 

    // Pass python values to cpp object.
    void update_c() {
        // Temporary containers
        Knots c_knots; // std::vector
        Coordinates c_control_points{}; // std::vector
        Degrees c_degrees{}; // std::array


        // Formulate Degrees
        py::buffer_info ds_buf = p_degrees.request();
        int* ds_buf_ptr = static_cast<int *>(ds_buf.ptr);
        for (i = 0; i < para_dim; i++) {
            c_degrees[i] = Degree{ds_buf_ptr[i]};
        }

        // Formulate Knotvectors
        i = 0;
        KnotVectors c_knot_vectors{};
        for (py::handle kv : p_knot_vectors) {
            c_knots.clear(); // Clear each time.
            for (py::handle k : kv) {
                c_knots.push_back(Knot(k.cast<double>()));
            }
            std::shared_ptr knot_vector{std::make_shared<KnotVector>(c_knots)};
            c_knot_vectors[i] = knot_vector;
            i++;
        }

        // Formulate Parameter_space
        c_parameter_space = std::make_shared<ParameterSpace>(c_knot_vectors, c_degrees);

        // Formulate Control Points
        c_control_points.clear();
        py::buffer_info cps_buf = p_control_points.request();
        double* cps_buf_ptr = static_cast<double *>(cps_buf.ptr);

        for (i = 0; i < cps_buf.shape[0]; i++) { // cps_buf.shape[0] : number of cps
            Coordinate control_point{};
            for (j = 0; j < dim; j++) { // dim : cps_buf.shape[1]
                control_point[j] = ScalarCoordinate{cps_buf_ptr[i * dim + j]};
            }
            c_control_points.push_back(control_point);
        }

        // Formulate Vector Space
        c_vector_space = std::make_shared<VectorSpace>(c_control_points);

        // Now, (re)initialize BSpline  
        c_bspline = BSpline{c_parameter_space, c_vector_space}; 
    }

    // Pass cpp object values to python.
    void update_p() {

        // Read from spline
        OutputInformation const &c_infos = c_bspline.Write();

        // Parameter space - knot vectors, degrees
        OutputParameterSpace const &parameter_space = std::get<0>(c_infos);
        OutputKnotVectors const &knot_vectors = std::get<0>(parameter_space);
        OutputDegrees const &degrees = std::get<1>(parameter_space);

        // Vector space - Coordinates(control points)
        OutputVectorSpace const &vector_space = std::get<1>(c_infos);
        OutputCoordinates const &coordinates = std::get<0>(vector_space);

        // Unpack - knot vectors
        p_knot_vectors.attr("clear")();
        for (auto& knotvector : knot_vectors) {
            py::list p_kv;
            for (auto& knot : knotvector) {
                p_kv.append(utilities::string_operations::ConvertToNumber<double>(knot));
            }
            p_knot_vectors.append(p_kv);
        }

        // Unpack - degrees
        // Edit existing array, since it never changes in size.
        p_degrees.resize({para_dim}); // inplace "flatten"
        py::buffer_info ds_buf = p_degrees.request();
        int* ds_buf_ptr = static_cast<int *>(ds_buf.ptr);

        i = 0;
        for (auto& degree : degrees) {
            ds_buf_ptr[i] = utilities::string_operations::ConvertToNumber<int>(degree);
            i++;
        }

        // Unpack - Coordinates (control points)
        p_control_points = py::array_t<double>(coordinates.size() * dim);
        py::buffer_info cps_buf = p_control_points.request();
        double* cps_buf_ptr = static_cast<double *>(cps_buf.ptr);

        i = 0;
        for (auto& coordinate : coordinates) {
            j = 0;
            for (auto& coord : coordinate) {
                cps_buf_ptr[i * dim + j] =
                    utilities::string_operations::ConvertToNumber<double>(coord);
                j++;
            }
            i++;
        }

        p_control_points.resize({(int) coordinates.size(), dim});

    }

    // Evaluate.
    py::array_t<double> evaluate(py::array_t<double> queries) {

        if (!skip_update) {
            update_c();
        }

        // Extract input array info.
        py::buffer_info q_buf = queries.request();
        double* q_buf_ptr = static_cast<double *>(q_buf.ptr);

        // Prepare results array.
        py::array_t<double> results(q_buf.shape[0] * dim);
        py::buffer_info r_buf = results.request();
        double* r_buf_ptr = static_cast<double *>(r_buf.ptr);

        // Loop.
        int num_queries = q_buf.shape[0];
        for (i = 0; i < num_queries; i++){
            ParametricCoordinate pc{};
            for (j = 0; j < para_dim; j++) {
                pc[j] = ScalarParametricCoordinate{(q_buf_ptr[i * para_dim + j])};
            }
            Coordinate const &c_result = c_bspline(pc);

            j = 0;
            for (auto& sc : c_result) { // `sc` : ScarlarCoordinate
                r_buf_ptr[i * dim + j] = sc.Get();
                j++;
            }
        }

        results.resize({num_queries, dim});

        return results;
    }

    // Derivative.
    py::array_t<double> derivative(py::array_t<double> queries, py::array_t<int> orders) {

        if (!skip_update) {
            update_c();
        }

        // Extract input arrays info.
        py::buffer_info q_buf = queries.request(), o_buf = orders.request();
        double* q_buf_ptr = static_cast<double *>(q_buf.ptr);
        int* o_buf_ptr = static_cast<int *>(o_buf.ptr);

        // Prepare results array.
        py::array_t<double> results(q_buf.shape[0] * dim);
        py::buffer_info r_buf = results.request();
        double* r_buf_ptr = static_cast<double *>(r_buf.ptr);

        // Formulate Derivative Orders.
        Derivative derivative{};
        for (i = 0; i < o_buf.shape[0]; i++) {
            derivative[i] = splinelib::Derivative{o_buf_ptr[i]};
        }

        // Loop - Queries.
        int num_queries = q_buf.shape[0];
        for (i = 0; i < num_queries; i++) {
            ParametricCoordinate pc{};
            for (j = 0; j < para_dim; j++) {
                pc[j] = ScalarParametricCoordinate{q_buf_ptr[i * para_dim + j]};
            }
            Coordinate const &c_result = c_bspline(pc, derivative);

            // Write `c_result` to `results`.
            j = 0;
            for (const auto& sc : c_result) { // `sc` : ScalarCoordinate
                r_buf_ptr[i * dim + j] = sc.Get();
                j++;
            }
        }

        results.resize({num_queries, dim});

        return results;

    }

    void insert_knots(int p_dim, py::list knots) {

        if (!skip_update) {
            update_c();
        }

        splinelib::Dimension inserting_p_dim{p_dim};
        for (py::handle k : knots) {
            c_bspline.InsertKnot(inserting_p_dim, Knot{k.cast<double>()});
        }
        update_p();
    }

    void remove_knots(int p_dim, py::list knots, double tol) {
        /*
         *  BSpline::RemoveKnot returns multiplicity. 
         *  In Debug mode, it returns -1 if you can't remove the knot or there's
         *    nothing to remove. At the same time in Debug mode, `splinelib::Multiplicity`
         *    is not supposed to be negative (see named_type.inc, line 25): you get an error.
         *  In Release mode, however, you get 1 even if you can't remove the knot
         *    or there's nothing to remove. Thus, this function will return nothing.
         *  Use with caution. You've been warned.
         */

        if (!skip_update) {
            update_c();
        }

        splinelib::Dimension removing_p_dim{p_dim};
        splines::Tolerance tolerance{tol};

        for (py::handle k : knots) {
            c_bspline.RemoveKnot(removing_p_dim, Knot{k.cast<double>()}, tolerance);
        }

        update_p();

    }
    
    void elevate_degree(int p_dim) {

        if (!skip_update) {
            update_c();
        }

        splinelib::Dimension elevating_p_dim{p_dim};
        c_bspline.ElevateDegree(elevating_p_dim);

        update_p();
    }

    bool reduce_degree(int p_dim, double tol) {

        if (!skip_update) {
            update_c();
        }

        bool reduced;

        splinelib::Dimension reducing_p_dim{p_dim};
        splines::Tolerance tolerance{tol};
        reduced = c_bspline.ReduceDegree(reducing_p_dim, tolerance);

        update_p();

        return reduced;
    }


    py::array_t<double> sample(py::array_t<int> query_resolutions) {

        if (!skip_update) {
            update_c();
        }

        // Extract input array info.
        py::buffer_info q_buf = query_resolutions.request();
        int* q_buf_ptr = static_cast<int *>(q_buf.ptr);

        // Prepare results array.
        int num_results = 1;
        for (i = 0; i < para_dim; i++) {
            num_results *= q_buf_ptr[i];
        }
        py::array_t<double> results(num_results * dim);
        py::buffer_info r_buf = results.request();
        double* r_buf_ptr = static_cast<double *>(r_buf.ptr);

        // Prepare NumberOfParametricCoordinates
        //   Could be done with "Prepare results array", but easier to read this way.
        NumberOfParametricCoordinates npc{};
        for (i = 0; i < para_dim; i ++) {
            npc[i] = splinelib::Length{q_buf_ptr[i]};
        }

        // Sample and write to `results`
        Coordinates sampled_coordinates = c_bspline.Sample(npc);
        for (int i = 0; i < sampled_coordinates.size(); i++) {
            Coordinate c = sampled_coordinates[i];

            j = 0;
            for (auto& sc : c) {
                r_buf_ptr[i * dim + j] = sc.Get();
                j++;
            }
        }

        results.resize({num_results, dim});

        return results;

    }

    double fit_curve(py::array_t<double> points,
                     int degree,
                     int num_control_points,
                     bool centripetal,
                     py::list knot_vectors) {

        if (para_dim != 1) {
            throw std::invalid_argument("Your spline's parametric dimension should be 1 for curve interpolation.");
        }

        // Extract input array info.
        py::buffer_info p_buf = points.request();
        double* p_buf_ptr = static_cast<double *>(p_buf.ptr);

        // Prepare vars for interpolation
        int num_points = p_buf.shape[0];
        int curve_dim = p_buf.shape[1];
        if (curve_dim != dim) {
            throw std::invalid_argument("Dimension mis-match between spline and interpolation query points.");
        }

        // Update KnotVector (is optional)
        std::vector<double> knot_vector, control_points;
        for (auto& knotvector : knot_vectors) {
            for (auto& knot : knotvector) {
                // Empty knot_vector is indicator for updating
                knot_vector.push_back(knot.cast<double>());
            }
        }

        double residual = FitCurve(p_buf_ptr, 
                                   num_points,
                                   curve_dim,
                                   degree,
                                   num_control_points,
                                   centripetal,
                                   knot_vector,
                                   control_points);

        // Write degree
        p_degrees = py::array_t<int>(para_dim);
        py::buffer_info pd_buf = p_degrees.request();
        int* pd_buf_ptr = static_cast<int *>(pd_buf.ptr);
        pd_buf_ptr[0] = degree;

        // Write knot vector
        p_knot_vectors.attr("clear")();
        py::list kv;
        for (const auto& k : knot_vector) {
            kv.append(k);
        }
        p_knot_vectors.append(kv);

        // Write control points
        p_control_points = py::array_t<double>(control_points.size());
        py::buffer_info pc_buf = p_control_points.request();
        double* pc_buf_ptr = static_cast<double *>(pc_buf.ptr);

        for (i = 0; i < num_control_points; i++) {
            for (j = 0; j < curve_dim; j++) {
                pc_buf_ptr[i * dim + j] = control_points[i * dim + j];
            }
        }

        p_control_points.resize({num_control_points, dim});

        update_c();

        return residual;

    }

    void interpolate_curve(py::array_t<double> points,
                           int degree,
                           bool centripetal) {

        py::buffer_info p_buf = points.request();
        int num_control_points = p_buf.shape[0];
        fit_curve(points, degree, num_control_points, centripetal, p_knot_vectors);

    }

    double approximate_curve(py::array_t<double> points,
                             int degree,
                             int num_control_points,
                             bool centripetal) {

        return fit_curve(points, degree, num_control_points, centripetal, p_knot_vectors);
    }

    void fit_surface(py::array_t<double> points,
                     int size_u,
                     int size_v,
                     int degree_u,
                     int degree_v,
                     bool centripetal) {

        if (para_dim != 2) {
            throw std::invalid_argument("Your spline's parametric dimension should be 2 for surface interpolation.");
        }

        // Extract input array info.
        py::buffer_info p_buf = points.request();
        double* p_buf_ptr = static_cast<double *>(p_buf.ptr);

        // Prepare vars for interpolation
        int num_points = p_buf.shape[0];
        int surface_dim = p_buf.shape[1];
        if (surface_dim != dim) {
            throw std::invalid_argument("Dimension mis-match between spline and interpolation query points.");
        }
        std::vector<double> knot_vector_u, knot_vector_v, control_points;

        FitSurface(p_buf_ptr,
                   num_points,
                   surface_dim,
                   degree_u,
                   degree_v,
                   size_u,
                   size_v,
                   centripetal,
                   knot_vector_u,
                   knot_vector_v,
                   control_points);

        // Write degree
        p_degrees = py::array_t<int>(para_dim);
        py::buffer_info pd_buf = p_degrees.request();
        int* pd_buf_ptr = static_cast<int *>(pd_buf.ptr);
        pd_buf_ptr[0] = degree_u;
        pd_buf_ptr[1] = degree_v;

        // Write knot vector
        p_knot_vectors.attr("clear")();
        py::list kv_u, kv_v;
        for (const auto& k : knot_vector_u) {
            kv_u.append(k);
        }
        for (const auto& k : knot_vector_v) {
            kv_v.append(k);
        }
        p_knot_vectors.append(kv_u);
        p_knot_vectors.append(kv_v);


        // Write control points
        p_control_points = py::array_t<double>(control_points.size());
        py::buffer_info pc_buf = p_control_points.request();
        double* pc_buf_ptr = static_cast<double *>(pc_buf.ptr);

        for (i = 0; i < num_points; i++) {
            for (j = 0; j < surface_dim; j++) {
                pc_buf_ptr[i * dim + j] = control_points[i * dim + j];
            }
        }

        p_control_points.resize({num_points, dim});

        update_c();

    }

    void interpolate_surface(py::array_t<double> points,
                             int size_u,
                             int size_v,
                             int degree_u,
                             int degree_v,
                             bool centripetal) {
        fit_surface(points, size_u, size_v, degree_u, degree_v, centripetal);
    }

    void write_iges(std::string fname) {

        input_output::iges::Write(
            {std::make_shared<BSpline>(c_parameter_space, c_vector_space)},
            fname
        );

    } 

    void write_xml(std::string fname) {

        input_output::xml::Write(
            {std::make_shared<BSpline>(c_parameter_space, c_vector_space)},
            fname
        );

    }

    void write_irit(std::string fname) {

        input_output::irit::Write(
            {std::make_shared<BSpline>(c_parameter_space, c_vector_space)},
            fname
        );

    }
};

template<int para_dim, int dim>
void add_bspline_pyclass(py::module &m, const char *class_name) {
    py::class_<PyBSpline<para_dim, dim>> klasse(m, class_name);

    klasse.def(py::init<>())
          .def(py::init<py::array_t<int>, py::list, py::array_t<double>>(),
                    py::arg("degrees"),
                    py::arg("knot_vectors"),
                    py::arg("control_points"))
          .def_readwrite("knot_vectors",
                             &PyBSpline<para_dim, dim>::p_knot_vectors)
          .def_readwrite("degrees",
                             &PyBSpline<para_dim, dim>::p_degrees)
          .def_readwrite("control_points",
                             &PyBSpline<para_dim, dim>::p_control_points)
          .def_readwrite("skip_update",
                             &PyBSpline<para_dim, dim>::skip_update)
          .def_readonly("whatami",
                             &PyBSpline<para_dim, dim>::whatami)
          .def("evaluate",
                   &PyBSpline<para_dim, dim>::evaluate,
                   py::arg("queries"),
                   py::return_value_policy::move)
          .def("derivative",
                   &PyBSpline<para_dim, dim>::derivative,
                   py::arg("queries"),
                   py::arg("orders"),
                   py::return_value_policy::move)
          .def("insert_knots",
                   &PyBSpline<para_dim, dim>::insert_knots,
                   py::arg("p_dim"),
                   py::arg("knots"))
          .def("remove_knots",
                   &PyBSpline<para_dim, dim>::remove_knots,
                   py::arg("p_dim"),
                   py::arg("knots"),
                   py::arg("tolerance"))
          .def("elevate_degree",
                   &PyBSpline<para_dim, dim>::elevate_degree,
                   py::arg("p_dim"))
          .def("reduce_degree",
                   &PyBSpline<para_dim, dim>::reduce_degree,
                   py::arg("p_dim"),
                   py::arg("tolerance"))
          .def("sample",
                   &PyBSpline<para_dim, dim>::sample,
                   py::arg("resoultion"),
                   py::return_value_policy::move)
          .def("write_iges",
                   &PyBSpline<para_dim, dim>::write_iges,
                   py::arg("fname"))
          .def("write_xml",
                   &PyBSpline<para_dim, dim>::write_xml,
                   py::arg("fname"))
          .def("write_irit",
                   &PyBSpline<para_dim, dim>::write_irit,
                   py::arg("fname"))
          .def("update_c",
                   &PyBSpline<para_dim, dim>::update_c)
          .def("update_p",
                   &PyBSpline<para_dim, dim>::update_p);


    if constexpr (para_dim == 1) {
        klasse.def("fit_curve",
                       &PyBSpline<para_dim, dim>::fit_curve,
                       py::arg("points"),
                       py::arg("degree"),
                       py::arg("num_control_points"),
                       py::arg("centripetal"),
                       py::arg("knot_vector"))
              .def("interpolate_curve",
                       &PyBSpline<para_dim, dim>::interpolate_curve,
                       py::arg("points"),
                       py::arg("degree"),
                       py::arg("centripetal"))
              .def("approximate_curve",
                       &PyBSpline<para_dim, dim>::approximate_curve,
                       py::arg("points"),
                       py::arg("degree"),
                       py::arg("num_control_points"),
                       py::arg("centripetal"));
    }

    if constexpr (para_dim == 2) {
        klasse.def("fit_surface",
                       &PyBSpline<para_dim, dim>::fit_surface,
                       py::arg("points"),
                       py::arg("size_u"),
                       py::arg("size_v"),
                       py::arg("degree_u"),
                       py::arg("degree_v"),
                       py::arg("centripetal"))
              .def("interpolate_surface",
                       &PyBSpline<para_dim, dim>::interpolate_surface,
                       py::arg("points"),
                       py::arg("size_u"),
                       py::arg("size_v"),
                       py::arg("degree_u"),
                       py::arg("degree_v"),
                       py::arg("centripetal"));
    }
}
