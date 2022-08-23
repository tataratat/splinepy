#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

// pybind11
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

// SplineLib
#include <Sources/InputOutput/iges.hpp>
#include <Sources/InputOutput/irit.hpp>
#include <Sources/InputOutput/vtk.hpp>
#include <Sources/InputOutput/xml.hpp>
#include <Sources/Splines/nurbs.hpp>
#include <Sources/Utilities/named_type.hpp>
#include <Sources/Utilities/string_operations.hpp>

// Local
#include <splinepy/py/py_rational_bezier.hpp>
#include <splinepy/splines/helpers.hpp>
#include <splinepy/splines/nurbs.hpp>

namespace py = pybind11;

using namespace splinelib::sources;

template <int para_dim, int dim>
class PyNurbs {
public:
  using Nurbs = splinepy::splines::Nurbs<para_dim, dim>;

  // For writing cpp splines
  using ParameterSpace = typename Nurbs::ParameterSpace_;
  using WeightedVectorSpace = typename Nurbs::WeightedVectorSpace_;
  using Coordinates = typename WeightedVectorSpace::Coordinates_;
  using Weights = typename WeightedVectorSpace::Weights_;
  using Weight = splinelib::Weight;
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
  using NumberOfParametricCoordinates =
      typename ParameterSpace::NumberOfParametricCoordinates_;

  int para_dim_ = para_dim;
  int dim_ = dim;

  // Counters
  int i = 0;
  int j = 0;

  // Fr. Nurbs herself
  Nurbs c_nurbs;
  std::shared_ptr<ParameterSpace> c_parameter_space;
  std::shared_ptr<WeightedVectorSpace> c_weighted_vector_space;

  // Fr. Nurbs' answer to the question "What am I?"
  const std::string whatami =
      "NURBS, parametric dimension: " + std::to_string(para_dim)
      + ", physical dimension: " + std::to_string(dim);

  // Fr. Nurbs' python Family
  py::list p_knot_vectors;
  py::array_t<int> p_degrees;
  py::array_t<double> p_control_points;
  py::array_t<double> p_weights;

  // Fr. Nurbs' switch to skip update_c
  bool skip_update = false;

  // Constructor.
  PyNurbs() {}
  PyNurbs(py::array_t<int> degrees,           //
          py::list knot_vectors,              //
          py::array_t<double> control_points, //
          py::array_t<double> weights)
      : p_degrees(degrees),
        p_knot_vectors(knot_vectors),
        p_control_points(control_points),
        p_weights(weights) {
    update_c();
  }

  // Pass python values to cpp object.
  void update_c() {
    // Temporary containers
    Knots c_knots;                  // std::vector
    Coordinates c_control_points{}; // std::vector
    Degrees c_degrees{};            // std::array
    Weights c_weights{};            // std::vector

    // Formulate Degrees
    py::buffer_info ds_buf = p_degrees.request();
    int* ds_buf_ptr = static_cast<int*>(ds_buf.ptr);
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
    c_parameter_space =
        std::make_shared<ParameterSpace>(c_knot_vectors, c_degrees);

    // Formulate Control Points
    py::buffer_info cps_buf = p_control_points.request();
    double* cps_buf_ptr = static_cast<double*>(cps_buf.ptr);

    c_control_points.clear();
    c_control_points.assign(cps_buf.shape[0], Coordinate{});

    for (i = 0; i < cps_buf.shape[0]; i++) { // cps_buf.shape[0] : number of cps
      Coordinate control_point{};
      for (j = 0; j < dim; j++) { // dim : cps_bus.shape[1]
        control_point[j] = ScalarCoordinate{cps_buf_ptr[i * dim + j]};
      }
      c_control_points[i] = control_point;
    }

    // Formulate Weights
    py::buffer_info ws_buf = p_weights.request();
    double* ws_buf_ptr = static_cast<double*>(ws_buf.ptr);

    c_weights.clear();
    c_weights.assign(ws_buf.shape[0], Weight{0.0});

    // Quick sanity check
    if (cps_buf.shape[0] != ws_buf.shape[0]) {
      throw std::invalid_argument(
          "Number of control points and number of weights does not match.");
    }

    for (i = 0; i < ws_buf.shape[0]; i++) {
      c_weights[i] = Weight{ws_buf_ptr[i]};
    }

    // Formulate Weighted Vector Space
    c_weighted_vector_space =
        std::make_shared<WeightedVectorSpace>(c_control_points, c_weights);

    // Now, (re)assign Fr. Nurbs
    c_nurbs = Nurbs{c_parameter_space, c_weighted_vector_space};
  }

  // Pass cpp object values to python.
  void update_p() {
    // ds
    p_degrees.resize({para_dim});
    py::buffer_info ds_buf = p_degrees.request();
    int* ds_buf_ptr = static_cast<int*>(ds_buf.ptr);
    c_nurbs.UpdateDegrees(ds_buf_ptr);

    // kvs
    c_nurbs.UpdateKnotVectors(p_knot_vectors);

    // cps & ws
    int numcps = c_nurbs.GetNCps();

    p_control_points = py::array_t<double>(numcps * dim);
    py::buffer_info cps_buf = p_control_points.request();
    double* cps_buf_ptr = static_cast<double*>(cps_buf.ptr);

    p_weights = py::array_t<double>(numcps);
    py::buffer_info ws_buf = p_weights.request();
    double* ws_buf_ptr = static_cast<double*>(ws_buf.ptr);

    c_nurbs.UpdateControlPointsAndWeights(cps_buf_ptr, ws_buf_ptr);

    p_control_points.resize({numcps, dim});
    p_weights.resize({numcps, 1});
  }

  // Evaluate.
  py::array_t<double> evaluate(py::array_t<double> queries) {
    if (!skip_update) {
      update_c();
    }

    // Extract input array info.
    py::buffer_info q_buf = queries.request();
    double* q_buf_ptr = static_cast<double*>(q_buf.ptr);

    // Prepare results array.
    py::array_t<double> results(q_buf.shape[0] * dim);
    py::buffer_info r_buf = results.request();
    double* r_buf_ptr = static_cast<double*>(r_buf.ptr);

    // Loop.
    int num_queries = q_buf.shape[0];
    for (i = 0; i < num_queries; i++) {
      ParametricCoordinate pc{};
      for (j = 0; j < para_dim; j++) {
        pc[j] = ScalarParametricCoordinate{(q_buf_ptr[i * para_dim + j])};
      }
      Coordinate const& c_result = c_nurbs(pc);
      j = 0;
      for (auto& sc : c_result) { // `sc` : ScarlarCoordinate
        r_buf_ptr[i * dim + j] = sc.Get();
        j++;
      }
    }

    results.resize({num_queries, dim});

    return results;
  }

  // multithread `evaluate` using std::thread
  py::array_t<double> p_evaluate(py::array_t<double> queries, int n_workers) {
    if (!skip_update) {
      update_c();
    }

    // Extract input array info.
    py::buffer_info q_buf = queries.request();
    double* q_buf_ptr = static_cast<double*>(q_buf.ptr);
    int n_queries = q_buf.shape[0];

    // Prepare results array.
    py::array_t<double> results(n_queries * dim);
    py::buffer_info r_buf = results.request();
    double* r_buf_ptr = static_cast<double*>(r_buf.ptr);

    // eval and fill up result array
    auto eval = [&](int begin, int end) {
      for (int id = begin; id < end; id++) { // n_queries
        ParametricCoordinate pc{};
        for (int pd = 0; pd < para_dim; pd++) {
          pc[pd] = ScalarParametricCoordinate{q_buf_ptr[id * para_dim + pd]};
        }
        Coordinate const& c_result = c_nurbs(pc);
        int d = 0;
        for (const auto& sc : c_result) {
          r_buf_ptr[id * dim + d] = sc.Get();
          d++;
        }
      }
    };

    // thread exe
    const int chunk_size = std::ceil(n_queries / n_workers);
    std::vector<std::thread> pool;
    for (i = 0; i < (n_workers - 1); i++) {
      std::thread th(eval, i * chunk_size, (i + 1) * chunk_size);
      pool.push_back(std::move(th));
    }
    {
      // last one
      std::thread th(eval, i * chunk_size, n_queries);
      pool.push_back(std::move(th));
    }

    for (auto& t : pool) {
      t.join();
    }

    results.resize({n_queries, dim});

    return results;
  }

  // Derivative.
  py::array_t<double> derivative(py::array_t<double> queries,
                                 py::array_t<int> orders) {
    if (!skip_update) {
      update_c();
    }

    // Extract input arrays info.
    py::buffer_info q_buf = queries.request(), o_buf = orders.request();
    double* q_buf_ptr = static_cast<double*>(q_buf.ptr);
    int* o_buf_ptr = static_cast<int*>(o_buf.ptr);

    // Prepare results array.
    py::array_t<double> results(q_buf.shape[0] * dim);
    py::buffer_info r_buf = results.request();
    double* r_buf_ptr = static_cast<double*>(r_buf.ptr);

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
      Coordinate const& c_result = c_nurbs(pc, derivative);

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

  // multithread `derivative` using std::thread
  py::array_t<double> p_derivative(py::array_t<double> queries, //
                                   py::array_t<int> orders,     //
                                   int n_workers) {
    if (!skip_update) {
      update_c();
    }

    // Extract input arrays info.
    py::buffer_info q_buf = queries.request(), o_buf = orders.request();
    double* q_buf_ptr = static_cast<double*>(q_buf.ptr);
    int* o_buf_ptr = static_cast<int*>(o_buf.ptr);
    int n_queries = q_buf.shape[0];

    // Prepare results array.
    py::array_t<double> results(n_queries * dim);
    py::buffer_info r_buf = results.request();
    double* r_buf_ptr = static_cast<double*>(r_buf.ptr);

    // Formulate Derivative Orders.
    Derivative derivative{};
    for (i = 0; i < o_buf.shape[0]; i++) {
      derivative[i] = splinelib::Derivative{o_buf_ptr[i]};
    }

    // deval and fill up result array
    auto deval = [&](int begin, int end) {
      for (int id = begin; id < end; id++) { // n_queries
        ParametricCoordinate pc{};
        for (int pd = 0; pd < para_dim; pd++) {
          pc[pd] = ScalarParametricCoordinate{q_buf_ptr[id * para_dim + pd]};
        }
        Coordinate const& c_result = c_nurbs(pc, derivative);
        int d = 0;
        for (const auto& sc : c_result) {
          r_buf_ptr[id * dim + d] = sc.Get();
          d++;
        }
      }
    };

    // thread exe
    const int chunk_size = std::ceil(n_queries / n_workers);
    std::vector<std::thread> pool;
    for (i = 0; i < (n_workers - 1); i++) {
      std::thread th(deval, i * chunk_size, (i + 1) * chunk_size);
      pool.push_back(std::move(th));
    }

    {
      // last one
      std::thread th(deval, i * chunk_size, n_queries);
      pool.push_back(std::move(th));
    }

    for (auto& t : pool) {
      t.join();
    }

    results.resize({n_queries, dim});

    return results;
  }

  /// Extract Elements to Rational Bezier Patches
  py::list ExtractBezierPatches() {
    const auto c_patches = splinepy::splines::ExtractBezierPatches(c_nurbs);
    py::list bezier_list{};
    for (std::size_t i{}; i < c_patches.size(); i++) {
      bezier_list.append(PyRationalBezier<para_dim, dim>{c_patches[i]});
    }
    return bezier_list;
  }

  // given parametric coordinate queires,
  // returns a tuple of (basis_functions, support_control_point_ids)
  py::tuple basis_functions(py::array_t<double> queries) {
    if (!skip_update) {
      update_c();
    }

    py::buffer_info q_buf = queries.request();
    double* q_buf_ptr = static_cast<double*>(q_buf.ptr);
    const int n_queries = q_buf.shape[0];

    // get prod(degrees. + 1)
    py::buffer_info ds_buf = p_degrees.request();
    int* ds_buf_ptr = static_cast<int*>(ds_buf.ptr);
    int n_supports = 1;
    for (i = 0; i < para_dim; i++) {
      n_supports *= (ds_buf_ptr[i] + 1);
    }

    // prepare results array
    py::array_t<double> basis_fn(n_queries * n_supports);
    py::buffer_info bf_buf = basis_fn.request();
    double* bf_buf_ptr = static_cast<double*>(bf_buf.ptr);

    py::array_t<int> support_cp_id(n_queries * n_supports);
    py::buffer_info sci_buf = support_cp_id.request();
    int* sci_buf_ptr = static_cast<int*>(sci_buf.ptr);

    // iterate queries and fill up result arrays
    for (i = 0; i < n_queries; i++) {
      ParametricCoordinate pc{};
      for (j = 0; j < para_dim; j++) {
        pc[j] = ScalarParametricCoordinate{q_buf_ptr[i * para_dim + j]};
      }
      double* bf_current_ptr = &bf_buf_ptr[i * n_supports];
      int* sci_current_ptr = &sci_buf_ptr[i * n_supports];
      c_nurbs.BasisFunctionsAndIDs(pc, bf_current_ptr, sci_current_ptr);
    }

    basis_fn.resize({n_queries, n_supports});
    support_cp_id.resize({n_queries, n_supports});

    return py::make_tuple(basis_fn, support_cp_id);
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
     * Nurbs::RemoveKnot returns multiplicity.
     *  In Debug mode, it returns -1 if you can't remove the knot or there's
     *    nothing to remove. At the same time in Debug mode,
     * `splinelib::Multiplicity` is not supposed to be negative (see
     * named_type.inc, line 25): you get an error. In Release mode, however, you
     * get 1 even if you can't remove the knot or there's nothing to remove.
     * Thus, this function will return nothing. Use with caution. You've been
     * warned.
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
    int* q_buf_ptr = static_cast<int*>(q_buf.ptr);

    // Prepare results array.
    int num_results = 1;
    for (i = 0; i < para_dim; i++) {
      num_results *= q_buf_ptr[i];
    }

    py::array_t<double> results(num_results * dim);
    py::buffer_info r_buf = results.request();
    double* r_buf_ptr = static_cast<double*>(r_buf.ptr);

    // Prepare NumberOfParametricCoordinates
    //   Could be done with "Prepare results array", but easier to read this
    //   way.
    NumberOfParametricCoordinates npc{};
    for (i = 0; i < para_dim; i++) {
      npc[i] = splinelib::Length{q_buf_ptr[i]};
    }

    // Sample and write to `results`
    Coordinates sampled_coordinates = c_nurbs.Sample(npc);
    for (int i = 0; i < sampled_coordinates.size(); i++) {
      Coordinate c = sampled_coordinates[i];

      j = 0;
      for (const auto& sc : c) {
        r_buf_ptr[i * dim + j] = sc.Get();
        j++;
      }
    }

    results.resize({num_results, dim});

    return results;
  }

  py::array_t<double> nearest_pcoord_midpoint(py::array_t<double> queries,
                                              int nthread) {
    // input arr
    py::buffer_info q_buf = queries.request();
    double* q_buf_ptr = static_cast<double*>(q_buf.ptr);
    const int n_queries = q_buf.shape[0];

    // output arr
    py::array_t<double> results(n_queries * para_dim);
    py::buffer_info r_buf = results.request();
    double* r_buf_ptr = static_cast<double*>(r_buf.ptr);

    // get proximity helper
    auto& proximity = c_nurbs.GetProximity();

    // lambda for nearest search
    auto nearest = [&](int begin, int end) {
      for (int id{begin}; id < end; id++) { // n_queries
        // ask the proximity guy
        const auto paracoord = proximity.FindNearestParametricCoordinate(
            &q_buf_ptr[id * dim],
            Nurbs::Proximity_::InitialGuess::MidPoint);
        // unpack
        for (int l{}; l < para_dim; ++l) {
          r_buf_ptr[id * para_dim + l] = paracoord[l].Get();
        }
      }
    };

    // nthreadexecution
    splinepy::utils::NThreadExecution(nearest, n_queries, nthread);

    results.resize({n_queries, para_dim});
    return results;
  }

  py::array_t<double> nearest_pcoord_kdt(py::array_t<double> queries,
                                         py::array_t<int> resolutions,
                                         int nthread) {
    // input arr
    py::buffer_info q_buf = queries.request();
    double* q_buf_ptr = static_cast<double*>(q_buf.ptr);
    const int n_queries = q_buf.shape[0];

    // resolution input
    int* qres_ptr = static_cast<int*>(resolutions.request().ptr);

    // output arr
    py::array_t<double> results(n_queries * para_dim);
    py::buffer_info r_buf = results.request();
    double* r_buf_ptr = static_cast<double*>(r_buf.ptr);

    // prepare proximity
    auto& proximity = c_nurbs.GetProximity();
    bool plant_newtree_please = true;
    // if resolutions of any entry is negative, we don't build a new tree
    std::array<int, para_dim> qres;
    for (int k{0}; k < para_dim; k++) {
      if (qres_ptr[k] < 0) {
        plant_newtree_please = false;
        break;
      }
      qres[k] = qres_ptr[k];
    }
    if (plant_newtree_please) {
      proximity.PlantNewKdTree(qres, nthread);
    }

    // lambda for nearest search
    auto nearest = [&](int begin, int end) {
      for (int id{begin}; id < end; id++) { // n_queries
        // ask the proximity guy
        const auto paracoord = proximity.FindNearestParametricCoordinate(
            &q_buf_ptr[id * dim],
            Nurbs::Proximity_::InitialGuess::KdTree);
        // unpack
        for (int l{}; l < para_dim; ++l) {
          r_buf_ptr[id * para_dim + l] = paracoord[l].Get();
        }
      }
    };

    // nthreadexecution
    splinepy::utils::NThreadExecution(nearest, n_queries, nthread);

    results.resize({n_queries, para_dim});
    return results;
  }

  void write_iges(std::string fname) {
    input_output::iges::Write(
        {std::make_shared<Nurbs>(c_parameter_space, c_weighted_vector_space)},
        fname);
  }

  void write_xml(std::string fname) {
    input_output::xml::Write(
        {std::make_shared<Nurbs>(c_parameter_space, c_weighted_vector_space)},
        fname);
  }

  void write_irit(std::string fname) {
    input_output::irit::Write(
        {std::make_shared<Nurbs>(c_parameter_space, c_weighted_vector_space)},
        fname);
  }
};

template <int para_dim, int dim>
void add_nurbs_pyclass(py::module& m, const char* class_name) {
  py::class_<PyNurbs<para_dim, dim>> klasse(m, class_name);

  klasse.def(py::init<>())
      .def(py::init<py::array_t<int>,
                    py::list,
                    py::array_t<double>,
                    py::array_t<double>>(),
           py::arg("degrees"),
           py::arg("knot_vectors"),
           py::arg("control_points"),
           py::arg("weights"))
      .def_readwrite("knot_vectors", &PyNurbs<para_dim, dim>::p_knot_vectors)
      .def_readwrite("degrees", &PyNurbs<para_dim, dim>::p_degrees)
      .def_readwrite("control_points",
                     &PyNurbs<para_dim, dim>::p_control_points)
      .def_readwrite("weights", &PyNurbs<para_dim, dim>::p_weights)
      .def_readwrite("skip_update", &PyNurbs<para_dim, dim>::skip_update)
      .def_readonly("whatami", &PyNurbs<para_dim, dim>::whatami)
      .def_readonly("dim", &PyNurbs<para_dim, dim>::dim_)
      .def_readonly("para_dim", &PyNurbs<para_dim, dim>::para_dim_)
      .def("evaluate",
           &PyNurbs<para_dim, dim>::evaluate,
           py::arg("queries"),
           py::return_value_policy::move)
      .def("derivative",
           &PyNurbs<para_dim, dim>::derivative,
           py::arg("queries"),
           py::arg("orders"),
           py::return_value_policy::move)
      .def("p_evaluate",
           &PyNurbs<para_dim, dim>::p_evaluate,
           py::arg("queries"),
           py::arg("n_threads"))
      .def("p_derivative",
           &PyNurbs<para_dim, dim>::p_derivative,
           py::arg("queries"),
           py::arg("orders"),
           py::arg("n_threads"))
      .def("basis_functions",
           &PyNurbs<para_dim, dim>::basis_functions,
           py::arg("queries"))
      .def("extract_bezier_patches",
           &PyNurbs<para_dim, dim>::ExtractBezierPatches)
      .def("insert_knots",
           &PyNurbs<para_dim, dim>::insert_knots,
           py::arg("p_dim"),
           py::arg("knots"))
      .def("remove_knots",
           &PyNurbs<para_dim, dim>::remove_knots,
           py::arg("p_dim"),
           py::arg("knots"),
           py::arg("tolerance"))
      .def("elevate_degree",
           &PyNurbs<para_dim, dim>::elevate_degree,
           py::arg("p_dim"))
      .def("reduce_degree",
           &PyNurbs<para_dim, dim>::reduce_degree,
           py::arg("p_dim"),
           py::arg("tolerance"))
      .def("sample",
           &PyNurbs<para_dim, dim>::sample,
           py::arg("resoultion"),
           py::return_value_policy::move)
      .def("nearest_pcoord_midpoint",
           &PyNurbs<para_dim, dim>::nearest_pcoord_midpoint,
           py::arg("queries"),
           py::arg("nthreads"))
      .def("nearest_pcoord_kdt",
           &PyNurbs<para_dim, dim>::nearest_pcoord_kdt,
           py::arg("queries"),
           py::arg("resolutions"),
           py::arg("nthreads"))
      .def("write_iges", &PyNurbs<para_dim, dim>::write_iges, py::arg("fname"))
      .def("write_xml", &PyNurbs<para_dim, dim>::write_xml, py::arg("fname"))
      .def("write_irit", &PyNurbs<para_dim, dim>::write_irit, py::arg("fname"))
      .def("update_c", &PyNurbs<para_dim, dim>::update_c)
      .def("update_p", &PyNurbs<para_dim, dim>::update_p)
      .def(py::pickle(
          [](const PyNurbs<para_dim, dim>& nurbs) {
            return py::make_tuple(nurbs.p_degrees,
                                  nurbs.p_knot_vectors,
                                  nurbs.p_control_points,
                                  nurbs.p_weights);
          },
          [](py::tuple t) {
            if (t.size() != 4) {
              throw std::runtime_error("Invalid PyNURBS state!");
            }

            PyNurbs<para_dim, dim> pyn(
                t[0].cast<py::array_t<int>>(),    // degrees
                t[1].cast<py::list>(),            // knot vectors
                t[2].cast<py::array_t<double>>(), // control points
                t[3].cast<py::array_t<double>>()  // weights
            );
            return pyn;
          }));
};
