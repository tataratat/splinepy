#pragma once

#include <string>

// pybind
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

//
#include "splinepy/py/py_spline.hpp"
#include "splinepy/utils/default_initialization_allocator.hpp"

namespace splinepy::py {

namespace py = pybind11;

// alias for frequently used vectors
using CoreSplineVector =
    std::vector<std::shared_ptr<splinepy::splines::SplinepyBase>>;
using IntVector = splinepy::utils::DefaultInitializationVector<int>;
using IntVectorVector = splinepy::utils::DefaultInitializationVector<IntVector>;
using DoubleVector = splinepy::utils::DefaultInitializationVector<double>;

/// @brief Extracts CoreSpline from list of PySplines.
/// @param pysplines
/// @param nthreads
/// @return
std::vector<PySpline::CoreSpline_> ToCoreSplineVector(py::list pysplines,
                                                      const int nthreads = 1);

/// @brief Overload to allow None entries. Will be filled with null spline
/// @param pysplines
/// @param para_dim_if_none
/// @param dim_if_none
/// @param nthreads
/// @return
std::vector<PySpline::CoreSpline_>
ToCoreSplineVector(py::list pysplines,
                   const int para_dim_if_none,
                   const int dim_if_none,
                   const int nthreads);

/// @brief Converts a CoreSplineVector to a Python spline list
/// @param splist
/// @return py::list
py::list ToPySplineList(CoreSplineVector& splist);

/*
 * Sort Vector using lambda expressions
 * https://stackoverflow.com/a/12399290
 */
template<typename T>
IntVector ArgSort(const std::vector<T>& v);

/// @brief raises if there's any mismatch between specified properties and all
/// the entries in the vector.
/// @param splist
/// @param name
/// @param para_dim
/// @param dim
/// @param degrees
/// @param control_mesh_resolutions
/// @param nthreads
void RaiseMismatch(const CoreSplineVector& splist,
                   const std::string name,
                   const int para_dim,
                   const int dim,
                   const IntVector& degrees,
                   const IntVector& control_mesh_resolutions,
                   const int nthreads);

/// @brief raises if there's any mismatch between two given vectors of splines.
/// @param splist0
/// @param splist1
/// @param name
/// @param para_dim
/// @param dim
/// @param degrees
/// @param control_mesh_resolutions
/// @param nthreads
void RaiseMismatch(const CoreSplineVector& splist0,
                   const CoreSplineVector& splist1,
                   const bool name,
                   const bool para_dim,
                   const bool dim,
                   const bool degrees,
                   const bool control_mesh_resolutions,
                   const int nthreads);

/**
 * @brief Determine the connectivity from center-vertices, assuming nothing of
 * the underlying grid
 *
 * Duplicate Points are not eliminated, assuming that a maximum of two points
 * are equivalent. If this is not the case an exception is thrown. In theory
 * this has complexity O(nlogn) whereas a KDTree has complexity O(n (logn)^dim).
 *
 * @tparam PhysicalPointType Type of Point coordinates
 * @tparam ScalarType Type determining the precision
 * @tparam parametric_dimension dimension of the object (e.g. surface in 3D)
 * @tparam boolean check_orientation to check if neighboring elements match
 *                          structured grid
 * @param face_center_vertices vertices in the centers of spline-surfaces
 * @param metric used for preordering the vertices along a line
 * @param tolerance tolerance (distance between two vertices that are joined)

 * @return connectivity as a std::vector<std::array<...>>
 */
py::array_t<int>
FindConnectivityFromCenters(const py::array_t<double>& face_center_vertices,
                            const int parametric_dimension,
                            const py::array_t<double>& metric,
                            const double tolerance);

/**
 * @brief  Determines the Connectivity of spline patches
 *
 * @param py_center_vertices  Vertices in the center of the boundaries
 * @param tolerance tolerance between two neighboring face centers for them
 * to be fused
 * @param parametric_dimension Parametric dimension of the spline grid
 * @return py::array_t<int> connectivity
 */
py::array_t<int>
InterfacesFromBoundaryCenters(const py::array_t<double>& py_center_vertices,
                              const double& tolerance,
                              const int& parametric_dimension);

/// @brief Orientation between two adjacent splines
///
/// If two splines share the same boundary this function retrieves their
/// orientation, by mapping the mappings of the parametric axis onto each other.
/// This is (among others) required for Gismo and Nutils export
///
/// @param pyspline_start Spline object from start
/// @param boundary_start Boundary ID from start spline
/// @param pyspline_end Spline object from end///to which is mapped
/// @param boundary_end Boundary ID of adjacent spline
/// @param tolerance Tolerance
/// @param int_mappings_ptr (output) integer mappings
/// @param bool_orientations_ptr (output) axis alignment
void GetBoundaryOrientation(
    const std::shared_ptr<splinepy::splines::SplinepyBase>& pyspline_start,
    const int& boundary_start,
    const std::shared_ptr<splinepy::splines::SplinepyBase>& pyspline_end,
    const int& boundary_end,
    const double& tolerance,
    int* int_mappings_ptr,
    bool* bool_orientations_ptr);

/// @brief Get the Boundary Orientations object
///
/// @param spline_list
/// @param base_id
/// @param base_face_id
/// @param neighbor_id
/// @param neighbor_face_id
/// @param tolerance
/// @param n_threads
/// @return py::tuple
py::tuple GetBoundaryOrientations(const py::list& spline_list,
                                  const py::array_t<int>& base_id,
                                  const py::array_t<int>& base_face_id,
                                  const py::array_t<int>& neighbor_id,
                                  const py::array_t<int>& neighbor_face_id,
                                  const double tolerance,
                                  const int n_threads);

/// @brief Extract all Boundary Patches and store them in a python list
///
/// @param spline_list List of splines
/// @param interfaces interfaces, with negative values for boundary elements
/// @param n_threads Number of threads
/// @return py::list
py::list ExtractAllBoundarySplines(const py::list& spline_list,
                                   const py::array_t<int>& interfaces,
                                   const int& n_threads);

/**
 * @brief  Adds a Boundary using a seed and G-continuity on boundary-splines
 *
 * This function might be a slight overkill, as it assigns all functions an ID,
 * even when previously assigned a different ID -> Future Project
 *
 * @param boundary_splines boundary patches
 * @param boundary_interfaces interfaces between boundary splines
 * @param global_interfaces global interfaces (in between "volume"-patches)
 * @param tolerance tolerance to be considered g1 (1 - cos(phi) < tolerance)
 * @param n_threads number of threads for parallel processing
 * @return int number of new boundaries
 */
int AddBoundariesFromContinuity(const py::list& boundary_splines,
                                const py::array_t<int>& boundary_interfaces,
                                py::array_t<int>& global_interfaces,
                                const double& tolerance,
                                const int& n_threads);

/// @brief Multi-patch splines. Here, use of the words
///        "patch" and "spline" are interchangeable
class PyMultipatch {
public:
  using CorePatches_ = std::vector<typename PySpline::CoreSpline_>;

  /// individual patches, a.k.a. splines
  py::list patches_;

  /// casted, core splines (patches) from patches_ (list of PySpline).
  CorePatches_ core_patches_;

  /// sub patches - similar concept to subelement / face elements
  std::shared_ptr<PyMultipatch> sub_multipatch_ = nullptr;

  /// Casted sub multi patch
  py::object py_sub_multipatch_ = py::none();

  /// center of each sub patch.
  py::array_t<double> sub_patch_centers_;

  /// ids for boundary_patches from sub_patches_
  py::array_t<int> boundary_patch_ids_;

  /// Boundary multi patches. Should consist of one smaller parametric dim
  std::shared_ptr<PyMultipatch> boundary_multipatch_ = nullptr;

  /// Casted boundary multi patches
  py::object py_boundary_multipatch_ = py::none();

  /// patch-to-patch connectivity information
  /// shape: (n_patches, n_boundary_elements)
  py::array_t<int> interfaces_;

  /// shape: (n_patches,)
  py::array_t<int> boundary_ids_;

  /// @brief Fields - they are saved as multi-patches
  py::list field_multipatches_;

  /// default number of threads for all the operations besides queries
  int n_default_threads_{1};

  /// hint flag that can be used to accelerate property decisions
  bool same_parametric_bounds_{false};

  /// internal flag to indicate that this a field multipatch that may contain
  /// null spline
  bool has_null_splines_{false};

  /// default tolerance
  double tolerance_{1e-11};

  /// ctor
  PyMultipatch() = default;

  /// move ctor
  PyMultipatch(PyMultipatch&& other) noexcept = default;

  /// @brief list (iterable) init -> pybind will cast to list for us if needed
  /// @param patches
  /// @param n_default_threads
  /// @param same_parametric_bounds
  PyMultipatch(py::list& patches,
               const int n_default_threads = 1,
               const bool same_parametric_bounds = false)
      : n_default_threads_(n_default_threads),
        same_parametric_bounds_(same_parametric_bounds) {
    SetPatchesNThreads(patches, n_default_threads);
  }

  /// @brief Internal use only. used to save field splines.
  /// @param patches
  /// @param n_default_threads
  /// @param para_dim_if_none
  /// @param dim_if_none
  PyMultipatch(py::list& patches,
               const int n_default_threads,
               const int para_dim_if_none,
               const int dim_if_none)
      : n_default_threads_(n_default_threads),
        has_null_splines_(true) {
    SetPatchesWithNullSplines(patches,
                              n_default_threads,
                              para_dim_if_none,
                              dim_if_none);
  }

  /// @brief given other shared pointer,  will copy references
  /// ctor. Expected use is for  to_derived call.
  ///
  /// @param other
  PyMultipatch(std::shared_ptr<PyMultipatch>& other) {
    patches_ = other->patches_;
    core_patches_ = other->core_patches_;
    sub_multipatch_ = other->sub_multipatch_;
    sub_patch_centers_ = other->sub_patch_centers_;
    py_sub_multipatch_ = other->py_sub_multipatch_;
    boundary_patch_ids_ = other->boundary_patch_ids_;
    boundary_multipatch_ = other->boundary_multipatch_;
    py_boundary_multipatch_ = other->py_boundary_multipatch_;
    interfaces_ = other->interfaces_;
    boundary_ids_ = other->boundary_ids_;
    field_multipatches_ = other->field_multipatches_;
    n_default_threads_ = other->n_default_threads_;
    same_parametric_bounds_ = other->same_parametric_bounds_;
    has_null_splines_ = other->has_null_splines_;
    tolerance_ = other->tolerance_;
  }

  /// @brief Gets core patches
  CorePatches_& CorePatches();

  /// @brief Clears
  void Clear();

  /// @brief returns para_dim of the first patch
  /// @return
  int ParaDim() { return CorePatches()[0]->SplinepyParaDim(); }

  /// @brief returns dim of the first patch
  /// @return
  int Dim() { return CorePatches()[0]->SplinepyDim(); }

  /// @brief my name is Multipatch
  /// @return
  std::string Name() { return "Multipatch"; }

  /// @brief Multipatch with n patches
  /// @return
  std::string WhatAmI() {
    return "Multipatch with " + std::to_string(core_patches_.size())
           + " patches.";
  }

  /// could overload, but it won't be clear for pybind
  /// @brief default patches setter, exposed to python
  /// @param patches
  void SetPatchesDefault(py::list& patches) {
    SetPatchesNThreads(patches, n_default_threads_);
  }

  /// @brief patches setter with nthreads
  /// @param patches
  /// @param nthreads
  void SetPatchesNThreads(py::list& patches, const int nthreads = 1);

  /// @brief patches setter for field patches that may contain null splines.
  /// @param patches
  /// @param nthreads
  /// @param para_dim_if_none
  /// @param dim_if_none
  void SetPatchesWithNullSplines(py::list& patches,
                                 const int nthreads,
                                 const int para_dim_if_none,
                                 const int dim_if_none);

  /// @brief Gets patches
  py::list GetPatches() { return patches_; }

  /// @brief returns sub patches as multi patches, uses default nthreads
  /// @return
  std::shared_ptr<PyMultipatch> SubMultipatch();

  /// @brief returns derived class of SubMultipatch
  /// @return
  py::object PySubMultipatch();

  /// @brief Computes centers of subpatches
  /// @return
  py::array_t<double> SubPatchCenters();

  /// @brief setter and getter for interface info. runs array shape check for
  /// setter
  /// @param interfaces
  py::array_t<int> Interfaces(const py::array_t<int>& interfaces);

  /// @brief Gets IDs of boundary patch
  py::array_t<int> BoundaryPatchIds();

  /// @brief Gets boundary multi patch
  std::shared_ptr<PyMultipatch> BoundaryMultipatch();

  /// @brief returns derived class of boundary multi patch
  /// @return
  py::object PyBoundaryMultipatch();

  /// @brief Evaluate at query points
  /// @param queries Query points
  /// @param nthreads Number of threads to use
  py::array_t<double> Evaluate(py::array_t<double> queries, const int nthreads);

  /// @brief Sample multi patch
  /// @param resolution
  /// @param nthreads Number of threads to use
  /// @param same_parametric_bounds
  py::array_t<double> Sample(const int resolution,
                             const int nthreads,
                             const bool same_parametric_bounds);

  /// @brief Adds fields
  /// @param fields
  /// @param check_name
  /// @param check_dims
  /// @param check_degrees
  /// @param check_control_mesh_resolutions
  /// @param nthreads
  void AddFields(py::list& fields,
                 const int field_dim,
                 const bool check_name,
                 const bool check_dims,
                 const bool check_degrees,
                 const bool check_control_mesh_resolutions,
                 const int nthreads);

  /// @brief Gets list of fields
  py::list GetFields() { return field_multipatches_; }

  /// @brief Get summed number of all control points
  int GetNumberOfControlPoints();

  /// @brief Get summed number of all control points
  py::array_t<int> GetControlPointOffsets();

  /// @brief Create a list of all control points
  py::array_t<double> GetControlPoints();
};

/// @brief ToDerived for PyMultipatches
/// @param core_obj
/// @return
py::object ToDerived(std::shared_ptr<PyMultipatch> core_obj);

} // namespace splinepy::py
