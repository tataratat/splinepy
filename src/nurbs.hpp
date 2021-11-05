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
#include <Sources/Splines/nurbs.hpp>
#include <Sources/Utilities/named_type.hpp>
#include <Sources/Utilities/string_operations.hpp>
#include <Sources/InputOutput/iges.hpp>
#include <Sources/InputOutput/xml.hpp>
#include <Sources/InputOutput/vtk.hpp>

namespace py = pybind11;

using namespace splinelib::sources;

template<int para_dim, int dim>
struct PyNurbs {

    using Nurbs = splines::Nurbs<para_dim, dim>;

    // For writing cpp splines
    using ParameterSpace = typename Nurbs::ParameterSpace_;
    using WeightedVectorSpace = typename Nurbs::WeightedVectorSpace_;
    using Coordinates = typename WeightedVectorSpace::Coordinates_;
    using Weights = typename WeightedVectorSpace::Weights_;
    using Weight = typename splinelib::Weight;
    using Degrees = typename ParameterSpace::Degrees_;
    using KnotVectors = typename ParameterSpace::KnotVectors_;
    using Coordinate = typename Coordinates::value_type;
    using Degree = typename Degrees::value_type;
    using KnotVector = typename KnotVectors::value_type::element_type;
    using Knots = typename KnotVector::Knots_;
    using ScalarCoordinate = typename Coordinate::value_type;
    using Knot = typename Knots::value_type;
    using KnotB = typename Nurbs::Knot_;
    using ParametricCoordinate = typename Nurbs::ParametricCoordinate_;
    using ScalarParametricCoordinate = typename ParametricCoordinate::value_type;
    using Derivative = typename Nurbs::Derivative_;
    using NumberOfParametricCoordinates = typename ParameterSpace::NumberOfParametricCoordinates_;

    // For reading cpp splines
    using OutputInformation = typename Nurbs::OutputInformation_;
    using OutputParameterSpace = typename std::tuple_element_t<0, OutputInformation>;
    using OutputWeightedVectorSpace = typename std::tuple_element_t<1, OutputInformation>;
    using OutputKnotVectors = typename std::tuple_element_t<0, OutputParameterSpace>;
    using OutputCoordinates = typename std::tuple_element_t<0, OutputWeightedVectorSpace>;
    using OutputWeights = typename std::tuple_element_t<1, OutputWeightedVectorSpace>;
    using OutputDegrees = typename std::tuple_element_t<1, OutputParameterSpace>;

    // For exporting cpp splines
    using ExportSplineItem = typename splines::SplineItem;
    using ExportSplineEntry = typename input_output::SplineEntry;
    using ExportSplines = typename input_output::Splines;

    // Counters
    int i = 0;
    int j = 0;

    // Fr. Nurbs herself
    Nurbs c_nurbs;
    std::shared_ptr<ParameterSpace> c_parameter_space;
    std::shared_ptr<WeightedVectorSpace> c_weighted_vector_space;

    // Fr. Nurbs' python Family
    py::list p_knot_vectors;
    py::array_t<int> p_degrees;
    py::array_t<double> p_control_points;
    py::array_t<double> p_weights;

    // Fr. Nurbs' switch to skip update_c
    bool skip_update = false;

    // Constructor.
    PyNurbs() {} 
    PyNurbs(py::array_t<int> degrees,
            py::array_t<double> control_points,
            py::list knot_vectors,
            py::array_t<double> weights) : 
                          p_degrees(degrees),
                          p_control_points(control_points),
                          p_knot_vectors(knot_vectors),
                          p_weights(weights) {
        update_c();
    } 

    // Pass python values to cpp object.
    void update_c() {
        // Temporary containers
        Knots c_knots; // std::vector
        Coordinates c_control_points{}; // std::vector
        Degrees c_degrees{}; // std::array
        Weights c_weights{}; // std::vector


        // Formulate Degrees
        std::fill(std::begin(c_degrees), std::end(c_degrees), Degree{0});
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
        py::buffer_info cps_buf = p_control_points.request();
        double* cps_buf_ptr = static_cast<double *>(cps_buf.ptr);

        c_control_points.clear();
        c_control_points.assign(cps_buf.shape[0], Coordinate{});

        for (i = 0; i < cps_buf.shape[0]; i++) { // cps_buf.shape[0] : number of cps
            Coordinate control_point{};
            for (j = 0; j < cps_buf.shape[1]; j++) { // cps_buf.shape[1] == dim 
                control_point[j] = ScalarCoordinate{cps_buf_ptr[i * dim + j]};
            }
            c_control_points[i] = control_point;
        }

        // Formulate Weights
        py::buffer_info ws_buf = p_weights.request();
        double* ws_buf_ptr = static_cast<double *>(ws_buf.ptr);

        c_weights.clear();
        c_weights.assign(ws_buf.shape[0], Weight{0.0});

        // Quick sanity check
        if (cps_buf.shape[0] != ws_buf.shape[0]) {
            throw std::invalid_argument(
                "Number of control points and number of weights does not match."
            );
        }

        for (i = 0; i < ws_buf.shape[0]; i++) {
            c_weights[i] = Weight{ws_buf_ptr[i]};
        }

        // Formulate Weighted Vector Space
        c_weighted_vector_space = std::make_shared<WeightedVectorSpace>(c_control_points, c_weights);

        // Now, (re)assign Fr. Nurbs 
        c_nurbs = Nurbs{c_parameter_space, c_weighted_vector_space}; 
    }

    // Pass cpp object values to python.
    void update_p() {

        // Read from spline
        OutputInformation const &c_infos = c_nurbs.Write();

        // Parameter space - knot vectors, degrees
        OutputParameterSpace const &parameter_space = std::get<0>(c_infos);
        OutputKnotVectors const &knot_vectors = std::get<0>(parameter_space);
        OutputDegrees const &degrees = std::get<1>(parameter_space);

        // Weighted Vector space - coordinates(control points), weights
        OutputWeightedVectorSpace const &weighted_vector_space = std::get<1>(c_infos);
        OutputCoordinates const &coordinates = std::get<0>(weighted_vector_space);
        OutputWeights const &weights = std::get<1>(weighted_vector_space);


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
        p_degrees = py::array_t<int>(para_dim);
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

        // Unpack - Weights 
        p_weights = py::array_t<double>(weights.size());
        py::buffer_info ws_buf = p_weights.request();
        double* ws_buf_ptr = static_cast<double *>(ws_buf.ptr);

        i = 0;
        for (auto& weight : weights) {
            ws_buf_ptr[i] =
                utilities::string_operations::ConvertToNumber<double>(weight);
            i++;
        }

        p_weights.resize({(int) weights.size(), 1}); // A tall vector

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
            Coordinate const &c_result = c_nurbs(pc);

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
            Coordinate const &c_result = c_nurbs(pc, derivative);

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
            c_nurbs.InsertKnot(inserting_p_dim, Knot{k.cast<double>()});
        }
        update_p();
    }

    void remove_knots(int p_dim, py::list knots, double tol) {
        /*
         *  Nurbs::RemoveKnot returns multiplicity. 
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
            c_nurbs.RemoveKnot(removing_p_dim, Knot{k.cast<double>()}, tolerance);
        }

        update_p();

    }
    
    void elevate_degree(int p_dim) {

        if (!skip_update) {
            update_c();
        }

        splinelib::Dimension elevating_p_dim{p_dim};
        c_nurbs.ElevateDegree(elevating_p_dim);

        update_p();
    }

    bool reduce_degree(int p_dim, double tol) {

        if (!skip_update) {
            update_c();
        }

        bool reduced;

        splinelib::Dimension reducing_p_dim{p_dim};
        splines::Tolerance tolerance{tol};
        reduced = c_nurbs.ReduceDegree(reducing_p_dim, tolerance);

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
        Coordinates sampled_coordinates = c_nurbs.Sample(npc);
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

    void write_iges(std::string fname) {

        input_output::iges::Write(
            {std::make_shared<Nurbs>(c_parameter_space, c_weighted_vector_space)},
            fname
        );

    } 

    void write_xml(std::string fname) {

        input_output::xml::Write(
            {std::make_shared<Nurbs>(c_parameter_space, c_weighted_vector_space)},
            fname
        );

    }

    void write_irit(std::string fname) {

        input_output::irit::Write(
            {std::make_shared<Nurbs>(c_parameter_space, c_weighted_vector_space)},
            fname
        );

    }
};
