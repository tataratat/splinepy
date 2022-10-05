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
#include <Sources/Splines/b_spline.hpp>
#include <Sources/Utilities/named_type.hpp>
#include <Sources/Utilities/string_operations.hpp>

// Local
#include <splinepy/fitting/fitting.hpp>
#include <splinepy/py/py_bezier.hpp>
#include <splinepy/splines/bspline.hpp>
#include <splinepy/splines/helpers.hpp>

namespace py = pybind11;

using namespace splinelib::sources;

template<int para_dim, int dim>
class PyBSpline {
public:
  using BSpline = splinepy::splines::BSpline<para_dim, dim>;

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
  using NumberOfParametricCoordinates =
      typename ParameterSpace::NumberOfParametricCoordinates_;

  // make both dims available
  int para_dim_ = para_dim;
  int dim_ = dim;

  // Hr. BSpline himself
  BSpline c_bspline;
  std::shared_ptr<ParameterSpace> c_parameter_space;
  std::shared_ptr<VectorSpace> c_vector_space;

  // Hr. BSpline's answer to the question "What am I?"
  const std::string whatami =
      "BSpline, parametric dimension: " + std::to_string(para_dim)
      + ", physical dimension: " + std::to_string(dim);

  // Hr. BSpline's python Family
  py::array_t<int> p_degrees;
  py::array_t<double> p_control_points;
  py::list p_knot_vectors;

  // Hr. BSpline's switch to skip_update_c
  bool skip_update = false;

  // Empty constructor
  PyBSpline() {
    py::list empty_list;
    p_knot_vectors.attr("clear")();
    p_knot_vectors.append(empty_list);
  }

  // Constructor
  PyBSpline(py::array_t<int> degrees,
            py::list knot_vectors,
            py::array_t<double> control_points)
      : p_degrees(degrees),
        p_control_points(control_points),
        p_knot_vectors(knot_vectors) {
    update_c();
  }

  // Pass python values to cpp object.
  void update_c() {
    // Temporary containers
    Knots c_knots;                  // std::vector
    Coordinates c_control_points{}; // std::vector
    Degrees c_degrees{};            // std::array

    // Formulate Degrees
    py::buffer_info ds_buf = p_degrees.request();
    int* ds_buf_ptr = static_cast<int*>(ds_buf.ptr);
    for (int i_para_dim = 0; i_para_dim < para_dim; i_para_dim++) {
      c_degrees[i_para_dim] = Degree{ds_buf_ptr[i_para_dim]};
    }

    // Formulate Knotvectors
    int i_kv = 0;
    KnotVectors c_knot_vectors{};
    for (py::handle kv : p_knot_vectors) {
      c_knots.clear(); // Clear each time.
      for (py::handle k : kv) {
        c_knots.push_back(Knot(k.cast<double>()));
      }
      std::shared_ptr knot_vector{std::make_shared<KnotVector>(c_knots)};
      c_knot_vectors[i_kv] = knot_vector;
      i_kv++;
    }

    // Formulate Parameter_space
    c_parameter_space =
        std::make_shared<ParameterSpace>(c_knot_vectors, c_degrees);

    // Formulate Control Points
    c_control_points.clear();
    py::buffer_info cps_buf = p_control_points.request();
    double* cps_buf_ptr = static_cast<double*>(cps_buf.ptr);

    for (int i_cps = 0; i_cps < cps_buf.shape[0];
         i_cps++) { // cps_buf.shape[0] : number of cps
      Coordinate control_point{};
      for (int j_dim = 0; j_dim < dim; j_dim++) { // dim : cps_buf.shape[1]
        control_point[j_dim] =
            ScalarCoordinate{cps_buf_ptr[i_cps * dim + j_dim]};
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
    // Unpack - knot vectors
    c_bspline.UpdateKnotVectors(p_knot_vectors);

    // Unpack - degrees
    // Edit existing array, since it never changes in size.
    p_degrees.resize({para_dim}); // inplace "flatten"
    py::buffer_info ds_buf = p_degrees.request();
    int* ds_buf_ptr = static_cast<int*>(ds_buf.ptr);
    c_bspline.UpdateDegrees(ds_buf_ptr);

    // Unpack - Coordinates (control points)
    int numcps = c_bspline.GetNCps();
    p_control_points = py::array_t<double>(numcps * dim);
    py::buffer_info cps_buf = p_control_points.request();
    double* cps_buf_ptr = static_cast<double*>(cps_buf.ptr);
    c_bspline.UpdateControlPoints(cps_buf_ptr);
    p_control_points.resize({numcps, dim});
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
    for (int i_query = 0; i_query < num_queries; i_query++) {
      ParametricCoordinate pc{};
      for (int j_para_dim = 0; j_para_dim < para_dim; j_para_dim++) {
        pc[j_para_dim] = ScalarParametricCoordinate{
            (q_buf_ptr[i_query * para_dim + j_para_dim])};
      }
      Coordinate const& c_result = c_bspline(pc);

      int j_dim = 0;
      for (auto& sc : c_result) { // `sc` : ScarlarCoordinate
        r_buf_ptr[i_query * dim + j_dim] = sc.Get();
        j_dim++;
      }
    }

    results.resize({num_queries, dim});

    return results;
  }

  // multithread `evaluate` using std::thread
  py::array_t<double> p_evaluate(py::array_t<double> queries, int n_workers) {
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
        Coordinate const& c_result = c_bspline(pc);
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
    for (int i_thread = 0; i_thread < (n_workers - 1); i_thread++) {
      std::thread th(eval, i_thread * chunk_size, (i_thread + 1) * chunk_size);
      pool.push_back(std::move(th));
    }
    {
      // last one
      std::thread th(eval, (n_workers - 1) * chunk_size, n_queries);
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
    for (int i_para_dim = 0; i_para_dim < o_buf.shape[0]; i_para_dim++) {
      derivative[i_para_dim] = splinelib::Derivative{o_buf_ptr[i_para_dim]};
    }

    // Loop - Queries.
    int num_queries = q_buf.shape[0];
    for (int i_query = 0; i_query < num_queries; i_query++) {
      ParametricCoordinate pc{};
      for (int j_para_dim = 0; j_para_dim < para_dim; j_para_dim++) {
        pc[j_para_dim] = ScalarParametricCoordinate{
            q_buf_ptr[i_query * para_dim + j_para_dim]};
      }
      Coordinate const& c_result = c_bspline(pc, derivative);

      // Write `c_result` to `results`.
      int j_dim = 0;
      for (const auto& sc : c_result) { // `sc` : ScalarCoordinate
        r_buf_ptr[i_query * dim + j_dim] = sc.Get();
        j_dim++;
      }
    }

    results.resize({num_queries, dim});

    return results;
  }

  // multithread `derivative` using std::thread
  py::array_t<double> p_derivative(py::array_t<double> queries, //
                                   py::array_t<int> orders,     //
                                   int n_workers) {
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
    for (int i_para_dim = 0; i_para_dim < o_buf.shape[0]; i_para_dim++) {
      derivative[i_para_dim] = splinelib::Derivative{o_buf_ptr[i_para_dim]};
    }

    // deval and fill up result array
    auto deval = [&](int begin, int end) {
      for (int id = begin; id < end; id++) { // n_queries
        ParametricCoordinate pc{};
        for (int pd = 0; pd < para_dim; pd++) {
          pc[pd] = ScalarParametricCoordinate{q_buf_ptr[id * para_dim + pd]};
        }
        Coordinate const& c_result = c_bspline(pc, derivative);
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
    for (int i_thread = 0; i_thread < (n_workers - 1); i_thread++) {
      std::thread th(deval, i_thread * chunk_size, (i_thread + 1) * chunk_size);
      pool.push_back(std::move(th));
    }
    {
      // last one
      std::thread th(deval, (n_workers - 1) * chunk_size, n_queries);
      pool.push_back(std::move(th));
    }

    for (auto& t : pool) {
      t.join();
    }

    results.resize({n_queries, dim});

    return results;
  }

  /// Extract Elements to polynomial Bezier patches
  py::list ExtractBezierPatches() {
    const auto c_patches = splinepy::splines::ExtractBezierPatches(c_bspline);
    py::list bezier_list{};
    for (std::size_t i_patch{}; i_patch < c_patches.size(); i_patch++) {
      bezier_list.append(PyBezier<para_dim, dim>{c_patches[i_patch]});
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
    for (int i_para_dim = 0; i_para_dim < para_dim; i_para_dim++) {
      n_supports *= (ds_buf_ptr[i_para_dim] + 1);
    }

    // prepare results array
    py::array_t<double> basis_fn(n_queries * n_supports);
    py::buffer_info bf_buf = basis_fn.request();
    double* bf_buf_ptr = static_cast<double*>(bf_buf.ptr);

    py::array_t<int> support_cp_id(n_queries * n_supports);
    py::buffer_info sci_buf = support_cp_id.request();
    int* sci_buf_ptr = static_cast<int*>(sci_buf.ptr);

    // iterate queries and fill up result arrays
    for (int i_query = 0; i_query < n_queries; i_query++) {
      ParametricCoordinate pc{};
      for (int j_para_dim = 0; j_para_dim < para_dim; j_para_dim++) {
        pc[j_para_dim] = ScalarParametricCoordinate{
            q_buf_ptr[i_query * para_dim + j_para_dim]};
      }
      double* bf_current_ptr = &bf_buf_ptr[i_query * n_supports];
      int* sci_current_ptr = &sci_buf_ptr[i_query * n_supports];
      c_bspline.BasisFunctionsAndIDs(pc, bf_current_ptr, sci_current_ptr);
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
      c_bspline.InsertKnot(inserting_p_dim, Knot{k.cast<double>()});
    }
    update_p();
  }

  void remove_knots(int p_dim, py::list knots, double tol) {
    /*
     *  BSpline::RemoveKnot returns multiplicity.
     *  In Debug mode, it returns -1 if you can't remove the knot or there's
     *    nothing to remove. At the same time in Debug mode,
     *    `splinelib::Multiplicity` is not supposed to be negative
     *    (see named_type.inc, line 25): you get an error.
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
    int* q_buf_ptr = static_cast<int*>(q_buf.ptr);

    // Prepare results array.
    int num_results = 1;
    for (int i_para_dim = 0; i_para_dim < para_dim; i_para_dim++) {
      num_results *= q_buf_ptr[i_para_dim];
    }
    py::array_t<double> results(num_results * dim);
    py::buffer_info r_buf = results.request();
    double* r_buf_ptr = static_cast<double*>(r_buf.ptr);

    // Prepare NumberOfParametricCoordinates
    //   Could be done with "Prepare results array", but easier to read this
    //   way.
    NumberOfParametricCoordinates npc{};
    for (int i_para_dim = 0; i_para_dim < para_dim; i_para_dim++) {
      npc[i_para_dim] = splinelib::Length{q_buf_ptr[i_para_dim]};
    }

    // Sample and write to `results`
    Coordinates sampled_coordinates = c_bspline.Sample(npc);
    for (std::size_t i_sample = 0; i_sample < sampled_coordinates.size();
         i_sample++) {
      Coordinate c = sampled_coordinates[i_sample];

      int j_dim = 0;
      for (const auto& sc : c) {
        r_buf_ptr[i_sample * dim + j_dim] = sc.Get();
        j_dim++;
      }
    }

    results.resize({num_results, dim});

    return results;
  }

  double fit_curve(py::array_t<double> points, //
                   int degree,                 //
                   int num_control_points,     //
                   bool centripetal,           //
                   py::list knot_vectors) {
    if (para_dim != 1) {
      throw std::invalid_argument(
          "parametric dimension should be 1 for curve interpolation.");
    }

    // Extract input array info.
    py::buffer_info p_buf = points.request();
    double* p_buf_ptr = static_cast<double*>(p_buf.ptr);

    // Prepare vars for interpolation
    int num_points = p_buf.shape[0];
    int curve_dim = p_buf.shape[1];
    if (curve_dim != dim) {
      throw std::invalid_argument(
          "Dimension mis-match between spline and interpolation query points.");
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
    int* pd_buf_ptr = static_cast<int*>(pd_buf.ptr);
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
    double* pc_buf_ptr = static_cast<double*>(pc_buf.ptr);

    for (int i_cps = 0; i_cps < num_control_points; i_cps++) {
      for (int j_dim = 0; j_dim < curve_dim; j_dim++) {
        pc_buf_ptr[i_cps * dim + j_dim] = control_points[i_cps * dim + j_dim];
      }
    }

    p_control_points.resize({num_control_points, dim});

    update_c();

    return residual;
  }

  void interpolate_curve(py::array_t<double> points, //
                         int degree,                 //
                         bool centripetal) {
    py::buffer_info p_buf = points.request();
    int num_control_points = p_buf.shape[0];
    fit_curve(points, degree, num_control_points, centripetal, p_knot_vectors);
  }

  double approximate_curve(py::array_t<double> points,
                           int degree,
                           int num_control_points,
                           bool centripetal) {
    return fit_curve(points,
                     degree,
                     num_control_points,
                     centripetal,
                     p_knot_vectors);
  }

  void fit_surface(py::array_t<double> points, //
                   int size_u,                 //
                   int size_v,                 //
                   int degree_u,               //
                   int degree_v,               //
                   bool centripetal) {
    if (para_dim != 2) {
      throw std::invalid_argument(
          "parametric dimension should be 2 for surface interpolation.");
    }

    // Extract input array info.
    py::buffer_info p_buf = points.request();
    double* p_buf_ptr = static_cast<double*>(p_buf.ptr);

    // Prepare vars for interpolation
    int num_points = p_buf.shape[0];
    int surface_dim = p_buf.shape[1];
    if (surface_dim != dim) {
      throw std::invalid_argument(
          "Dimension mis-match between spline and interpolation query points.");
    }
    std::vector<double> knot_vector_u, knot_vector_v, control_points;

    FitSurface(p_buf_ptr,     //
               num_points,    //
               surface_dim,   //
               degree_u,      //
               degree_v,      //
               size_u,        //
               size_v,        //
               centripetal,   //
               knot_vector_u, //
               knot_vector_v, //
               control_points //
    );

    // Write degree
    p_degrees = py::array_t<int>(para_dim);
    py::buffer_info pd_buf = p_degrees.request();
    int* pd_buf_ptr = static_cast<int*>(pd_buf.ptr);
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
    double* pc_buf_ptr = static_cast<double*>(pc_buf.ptr);

    for (int i_cps = 0; i_cps < num_points; i_cps++) {
      for (int j_dim = 0; j_dim < surface_dim; j_dim++) {
        pc_buf_ptr[i_cps * dim + j_dim] = control_points[i_cps * dim + j_dim];
      }
    }

    p_control_points.resize({num_points, dim});

    update_c();
  }

  void interpolate_surface(py::array_t<double> points, //
                           int size_u,                 //
                           int size_v,                 //
                           int degree_u,               //
                           int degree_v,               //
                           bool centripetal) {
    fit_surface(points, size_u, size_v, degree_u, degree_v, centripetal);
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
    auto& proximity = c_bspline.GetProximity();

    // lambda for nearest search
    auto nearest = [&](int begin, int end) {
      for (int id{begin}; id < end; id++) { // n_queries
        // ask the proximity guy
        const auto paracoord = proximity.FindNearestParametricCoordinate(
            &q_buf_ptr[id * dim],
            BSpline::Proximity_::InitialGuess::MidPoint);
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
    auto& proximity = c_bspline.GetProximity();
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
            BSpline::Proximity_::InitialGuess::KdTree);
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
        {std::make_shared<BSpline>(c_parameter_space, c_vector_space)},
        fname);
  }

  void write_xml(std::string fname) {
    input_output::xml::Write(
        {std::make_shared<BSpline>(c_parameter_space, c_vector_space)},
        fname);
  }

  void write_irit(std::string fname) {
    input_output::irit::Write(
        {std::make_shared<BSpline>(c_parameter_space, c_vector_space)},
        fname);
  }
};

template<int para_dim, int dim>
void add_bspline_pyclass(py::module& m, const char* class_name) {
  py::class_<PyBSpline<para_dim, dim>> klasse(m, class_name);

  klasse.def(py::init<>())
      .def(py::init<py::array_t<int>, py::list, py::array_t<double>>(),
           py::arg("degrees"),
           py::arg("knot_vectors"),
           py::arg("control_points"))
      .def_readwrite("knot_vectors", &PyBSpline<para_dim, dim>::p_knot_vectors)
      .def_readwrite("degrees", &PyBSpline<para_dim, dim>::p_degrees)
      .def_readwrite("control_points",
                     &PyBSpline<para_dim, dim>::p_control_points)
      .def_readwrite("skip_update", &PyBSpline<para_dim, dim>::skip_update)
      .def_readonly("whatami", &PyBSpline<para_dim, dim>::whatami)
      .def_readonly("dim", &PyBSpline<para_dim, dim>::dim_)
      .def_readonly("para_dim", &PyBSpline<para_dim, dim>::para_dim_)
      .def("evaluate",
           &PyBSpline<para_dim, dim>::evaluate,
           py::arg("queries"),
           py::return_value_policy::move)
      .def("derivative",
           &PyBSpline<para_dim, dim>::derivative,
           py::arg("queries"),
           py::arg("orders"),
           py::return_value_policy::move)
      .def("p_evaluate",
           &PyBSpline<para_dim, dim>::p_evaluate,
           py::arg("queries"),
           py::arg("n_threads"))
      .def("p_derivative",
           &PyBSpline<para_dim, dim>::p_derivative,
           py::arg("queries"),
           py::arg("orders"),
           py::arg("n_threads"))
      .def("extract_bezier_patches",
           &PyBSpline<para_dim, dim>::ExtractBezierPatches)
      .def("basis_functions",
           &PyBSpline<para_dim, dim>::basis_functions,
           py::arg("queries"))
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
      .def("nearest_pcoord_midpoint",
           &PyBSpline<para_dim, dim>::nearest_pcoord_midpoint,
           py::arg("queries"),
           py::arg("nthreads"))
      .def("nearest_pcoord_kdt",
           &PyBSpline<para_dim, dim>::nearest_pcoord_kdt,
           py::arg("queries"),
           py::arg("resolutions"),
           py::arg("nthreads"))
      .def("write_iges",
           &PyBSpline<para_dim, dim>::write_iges,
           py::arg("fname"))
      .def("write_xml", &PyBSpline<para_dim, dim>::write_xml, py::arg("fname"))
      .def("write_irit",
           &PyBSpline<para_dim, dim>::write_irit,
           py::arg("fname"))
      .def("update_c", &PyBSpline<para_dim, dim>::update_c)
      .def("update_p", &PyBSpline<para_dim, dim>::update_p)
      .def(py::pickle(
          [](const PyBSpline<para_dim, dim>& bspl) {
            return py::make_tuple(bspl.p_degrees,
                                  bspl.p_knot_vectors,
                                  bspl.p_control_points);
          },
          [](py::tuple t) {
            if (t.size() != 3) {
              throw std::runtime_error("Invalid PyBspline state!");
            }

            PyBSpline<para_dim, dim> pyb(t[0].cast<py::array_t<int>>(),
                                         t[1].cast<py::list>(),
                                         t[2].cast<py::array_t<double>>());

            return pyb;
          }));

  if constexpr (para_dim == 1) {
    klasse
        .def("fit_curve",
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
    klasse
        .def("fit_surface",
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
