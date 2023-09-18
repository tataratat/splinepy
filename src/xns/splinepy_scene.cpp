#include <array>
#include <memory>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "splinepy/py/py_spline.hpp"
#include "splinepy/utils/default_initialization_allocator.hpp"
#include "splinepy/utils/nthreads.hpp"
#include "splinepy/utils/print.hpp"

namespace splinepy::xns {

namespace py = pybind11;

template<typename T>
using DVector = splinepy::utils::DefaultInitializationVector<T>;
using PySpline = splinepy::py::PySpline;

static int DIM{};
static int BOUNDARY_PARA_DIM{};
static std::shared_ptr<PySpline> SPLINE{};
static std::shared_ptr<PySpline> BOUNDARY_SPLINE{};
static std::array<double, 3> LOWER_BOUND;
static std::array<double, 3> UPPER_BOUND;

/// @brief array memory. intended for fast (temporary) use with
/// contiguous memory layout.
/// for 2d-like access, set `stride_` and `is_2d` true
///
/// for non-owing data, use default init and set data and size yourself.
/// size is only used to check bounds in debug mode.
///
/// @tparam Type
template<typename Type, int dim = -1>
struct Data {
  Type* data_{nullptr};
  int size_{0};
  int stride0_{1}; // set directly. in case of 3d this should be dim0 * dim1
  int stride1_{1}; // set directly.

  /// @brief can't set size twice
  /// @param n
  void SetSize(int const& n) {
    if (data_) {
      splinepy::utils::PrintAndThrowError(
          "data is not empty and can't assign new size.");
    }
    data_ = new Type[n];
    size_ = n;
  }

  Data() = default;
  Data(const int n) { SetSize(n); }

  /// n is size. n * stride is size * stride, don't be confused.
  Data(int n, const int stride) : Data(n), stride0_(stride) {
    assert(n % stride == 0);
  }
  ~Data() { delete[] data_; }

  constexpr Type* data() {
    assert(data_);

    return data_;
  }
  constexpr const Type* data() const {
    assert(data_);

    return data_;
  }
  constexpr int size() const { return size_; }

  constexpr Type* begin() {
    assert(data_);

    return data_;
  }

  constexpr Type* end() {
    assert(data_);

    return data_ + size_;
  }

  constexpr void fill(Type const& value) {
    std::fill(data_, data_ + size_, value);
  }

  template<typename IndexType>
  constexpr Type& operator[](const IndexType& i) {
    assert(i < size_);
    assert(data_);

    return data_[i];
  }
  template<typename IndexType>
  constexpr const Type& operator[](const IndexType& i) const {
    assert(i < size_);
    assert(data_);

    return data_[i];
  }
  template<typename IndexType>
  constexpr Type& operator()(const IndexType& i, const IndexType& j) {
    static_assert(dim == 2);
    assert(i * stride0_ + j < size_);
    assert(data_);

    return data_[i * stride0_ + j];
  }
  template<typename IndexType>
  constexpr const Type& operator()(const IndexType& i,
                                   const IndexType& j) const {
    static_assert(dim == 2);
    assert(i * stride0_ + j < size_);
    assert(data_);

    return data_[i * stride0_ + j];
  }
  template<typename IndexType>
  constexpr Type&
  operator()(const IndexType& i, const IndexType& j, const IndexType& k) {
    static_assert(dim == 3);
    assert(i * stride0_ + j * stride1_ + k < size_);
    assert(data_);

    return data_[i * stride0_ + j * stride1_ + k];
  }
  template<typename IndexType>
  constexpr const Type&
  operator()(const IndexType& i, const IndexType& j, const IndexType& k) const {
    static_assert(dim == 3);
    assert(i * stride0_ + j * stride1_ + k < size_);
    assert(data_);

    return data_[i * stride0_ + j * stride1_ + k];
  }

  template<typename IndexType>
  constexpr Type* Pointer(const IndexType& i) {
    assert(i < size_);
    assert(data_);

    return &data_[i];
  }
  template<typename IndexType>
  constexpr const Type* Pointer(const IndexType& i) const {
    assert(i < size_);
    assert(data_);

    return &data_[i];
  }
  template<typename IndexType>
  constexpr Type* Pointer(const IndexType& i, const IndexType& j) {
    static_assert(dim == 2);
    assert(i * stride0_ + j < size_);
    assert(data_);

    return &data_[i * stride0_ + j];
  }
  template<typename IndexType>
  constexpr const Type* Pointer(const IndexType& i, const IndexType& j) const {
    static_assert(dim == 2);
    assert(i * stride0_ + j < size_);
    assert(data_);

    return &data_[i * stride0_ + j];
  }
  template<typename IndexType>
  constexpr Type*
  Pointer(const IndexType& i, const IndexType& j, const IndexType& k) const {
    static_assert(dim == 3);
    assert(i * stride0_ + j * stride1_ + k < size_);
    assert(data_);

    return &data_[i * stride0_ + j * stride1_ + k];
  }
  template<typename IndexType>
  constexpr Type*
  Pointer(const IndexType& i, const IndexType& j, const IndexType& k) {
    static_assert(dim == 3);
    assert(i * stride0_ + j * stride1_ + k < size_);
    assert(data_);

    return &data_[i * stride0_ + j * stride1_ + k];
  }
};

/// normal - we define a rule here, and it'd be your job to prepare
/// the foreign spline in a correct orientation
template<int dim, bool unit_normal = true, typename ArrayType>
inline static void Normal(const Data<double, 2>& first_dir, ArrayType& normal) {
  assert(first_dir.size() == 2 || first_dir.size() == 6);

  if constexpr (dim == 2) {
    const double& d0 = first_dir[0];
    const double& d1 = first_dir[1];

    if constexpr (unit_normal) {
      const double inv_norm2 = 1. / std::sqrt(d0 * d0 + d1 * d1);

      normal[0] = d1 * inv_norm2;
      normal[1] = -d0 * inv_norm2;
    } else {
      normal[0] = d1;
      normal[1] = -d0;
    }

  } else if constexpr (dim == 3) {
    const double& d0 = first_dir[0];
    const double& d1 = first_dir[1];
    const double& d2 = first_dir[2];
    const double& d3 = first_dir[3];
    const double& d4 = first_dir[4];
    const double& d5 = first_dir[5];

    if constexpr (unit_normal) {
      const double n0 = d1 * d5 - d2 * d4;
      const double n1 = d2 * d3 - d0 * d5;
      const double n2 = d0 * d4 - d1 * d3;

      const double inv_norm2 = 1. / std::sqrt(n0 * n0 + n1 * n1 + n2 * n2);

      normal[0] = n0 * inv_norm2;
      normal[1] = n1 * inv_norm2;
      normal[2] = n2 * inv_norm2;

    } else {

      normal[0] = d1 * d5 - d2 * d4;
      normal[1] = d2 * d3 - d0 * d5;
      normal[2] = d0 * d4 - d1 * d3;
    }
  } else {
    static_assert(dim == 2 || dim == 3, "unsupported dim");
  }
}

inline static bool IsNormalGapPositive(const Data<double, 2>& first_der,
                                       const Data<double>& nearest_minus_query,
                                       Data<double>& normal) {

  double normal_distance;
  const int dim = normal.size();

  // let's get normal distance
  if (dim == 2) {
    Normal<2, true>(first_der, normal); // I think this doesn't need to be unit
    normal_distance =
        normal[0] * nearest_minus_query[0] + normal[1] * nearest_minus_query[1];
  } else {
    Normal<3, true>(first_der, normal);
    normal_distance = normal[0] * nearest_minus_query[0]
                      + normal[1] * nearest_minus_query[1]
                      + normal[2] * nearest_minus_query[2];
  }

  // by definition, normal gap is negative of normal distance
  return normal_distance < 0.; // same as -normal_distance > 0.;
}

static void SetAABB(py::array_t<double>& aabb) {
  splinepy::utils::PrintInfo("setting AABB");

  if (static_cast<int>(aabb.shape(1)) != DIM) {
    splinepy::utils::PrintAndThrowError("dim mismatch between aabb and dim");
  }

  const double* ptr = static_cast<double*>(aabb.request().ptr);
  for (int i{}; i < DIM; ++i) {
    LOWER_BOUND[i] = ptr[i];
    UPPER_BOUND[i] = ptr[i + DIM];
  }
}

static void SetSpline(const std::shared_ptr<PySpline>& spline) {
  splinepy::utils::PrintInfo("setting spline");
  // set spline
  SPLINE = spline;
  // set dim
  DIM = SPLINE->Core()->SplinepyDim();
  // find aabb
  // typename splinepy::splines::SplinepyBase::ControlPointPointers_ cp_ptrs;
  // if (SPLINE->Core()->SplinepyIsRational()) {
  //  cp_ptrs = SPLINE->Core()->SplinepyWControlPointPointers
  //} else {

  //}
}

static void SetBoundary(const std::shared_ptr<PySpline>& spline) {
  splinepy::utils::PrintInfo("setting bdr spline");
  BOUNDARY_SPLINE = spline;
  BOUNDARY_PARA_DIM = BOUNDARY_SPLINE->Core()->SplinepyParaDim();
}

static void IsInSpline(const double* queries,
                       int* true_false,
                       const int& n_queries,
                       const int& nthreads) {
  splinepy::utils::PrintInfo("is in spline");

  DVector<DVector<int>> in_aabb_per_thread(nthreads);
  Data<int> in_aabb(n_queries);

  auto check_aabb_contains =
      [&](const int begin, const int end, const int i_thread) {
        auto& t_in_aabb = in_aabb_per_thread[i_thread];
        t_in_aabb.reserve(end - begin);

        for (int i{begin}; i < end; ++i) {
          // let's set all the query to 1 first - they are all in per default
          true_false[i] = 1;
          const double* i_query = &queries[i * DIM];

          // check AABB
          bool inside_aabb{true};
          for (int j{}; j < DIM; ++j) {
            const auto& j_q = i_query[j];

            // is this outside?
            if (j_q < LOWER_BOUND[i] || j_q > UPPER_BOUND[i]) {
              // if any of the point is outside we can break the loop
              inside_aabb = false;
              true_false[i] = 0; // mark return array already
              break;             // break j loop
            }
          }
          // we mark this
          if (inside_aabb) {
            t_in_aabb.push_back(i);
          }
        }
      };

  splinepy::utils::NThreadExe(check_aabb_contains, n_queries, nthreads);

  // put them all together
  int in_aabb_total{};
  for (auto& iapt : in_aabb_per_thread) {
    const int this_size = static_cast<int>(iapt.size());
    std::copy_n(iapt.begin(), this_size, in_aabb.Pointer(in_aabb_total));
    in_aabb_total += this_size;
  }

  // make sure trees are planted before you call this one
  //
  // 2 options:
  //   - 1. call spline's proximity and check if the nearest point is inside the
  //        spline. If it is on the boundary, it is not in side the spline
  //   - 2. call boundary spline's proximity and check normal gap and distance
  // we will go for 2, it's cheaper
  auto check_normal_gap = [&](const int begin,
                              const int end,
                              const int i_thread) {
    Data<double> parametric(BOUNDARY_PARA_DIM);
    Data<double> physical(DIM);
    Data<double> physical_minus_query(DIM);
    double distance{};
    double convergence{};
    double tolerance{1e-12};
    int max_iter{20 * BOUNDARY_PARA_DIM};
    Data<double, 2> first_derivatives(BOUNDARY_PARA_DIM * DIM);
    Data<double> normal(DIM); // tmp array for n Gap.
    for (int i{begin}; i < end; ++i) {
      const int& query_id = in_aabb[i];
      const double* query = &queries[query_id * DIM];

      // query
      BOUNDARY_SPLINE->Core()->SplinepyVRDMUMQuery(query,
                                                   tolerance,
                                                   max_iter,
                                                   false,
                                                   parametric.data(),
                                                   physical.data(),
                                                   physical_minus_query.data(),
                                                   distance,
                                                   convergence,
                                                   first_derivatives.data());

      // if distance is zero, this is not in
      // if normal gap is positive, it is not in.
      if (distance < tolerance
          || IsNormalGapPositive(first_derivatives,
                                 physical_minus_query,
                                 normal)) {
        true_false[query_id] = 0.;
        continue;
      }
    }
  };

  splinepy::utils::NThreadExe(check_normal_gap, in_aabb_total, nthreads);
}

static void NearestBoundaryPoint(const double* queries,
                                 const int& n_queries,
                                 const int& nthreads,
                                 double* results) {
  splinepy::utils::PrintInfo("nearest boundary point");

  auto nearest_query = [&](const int begin, const int end, const int i_thread) {
    Data<double> parametric(BOUNDARY_PARA_DIM);
    Data<double> physical_minus_query(DIM);
    double distance{};
    double convergence{};
    double tolerance{1e-12};
    int max_iter{20 * BOUNDARY_PARA_DIM};
    Data<double, 2> first_derivatives(BOUNDARY_PARA_DIM * DIM); // umsonst

    for (int i{begin}; i < end; ++i) {

      // query - we really just need physical
      BOUNDARY_SPLINE->Core()->SplinepyVRDMUMQuery(&queries[i * DIM],
                                                   tolerance,
                                                   max_iter,
                                                   false,
                                                   parametric.data(),
                                                   &results[i * DIM],
                                                   physical_minus_query.data(),
                                                   distance,
                                                   convergence,
                                                   first_derivatives.data());
    }
  };

  splinepy::utils::NThreadExe(nearest_query, n_queries, nthreads);
}

} // namespace splinepy::xns

// export for fortran
extern "C" {

void is_in_spline(double* queries,
                  int* true_false,
                  int* n_queries,
                  int* nthreads) {
  splinepy::xns::IsInSpline(queries, true_false, *n_queries, *nthreads);
};

void nearest_boundary_point(double* queries,
                            int* n_queries,
                            int* nthreads,
                            double* results) {
  splinepy::xns::NearestBoundaryPoint(queries, *n_queries, *nthreads, results);
};
}
