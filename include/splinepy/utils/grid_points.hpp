#pragma once

#include <cmath>
#include <numeric>
#include <utility>

#include <splinepy/utils/default_initialization_allocator.hpp>
#include <splinepy/utils/print.hpp>

namespace splinepy::utils {

/// @brief GridPoints based dynamic grid point sampler
class GridPoints {
public:
  template<typename T>
  using Vector_ = splinepy::utils::DefaultInitializationVector<T>;

  /// @brief Vector of bounds
  Vector_<double> bounds_;
  /// @brief Vector of resolutions
  Vector_<int> resolutions_;
  /// @brief Vector of vectors of entries
  Vector_<Vector_<double>> entries_;
  /// @brief each dimension's interval
  Vector_<double> step_size_;
  /// @brief Length
  int len_;
  /// @brief Dimension
  int dim_;

  GridPoints() = default;
  GridPoints(const int dim, const double* bounds, const int* resolutions) {
    SetUp(dim, bounds, resolutions);
  }

  /// @brief Set up resolutions and entries
  /// @param dim
  /// @param bounds
  /// @param resolutions
  void SetUp(const int dim, const double* bounds, const int* resolutions) {
    dim_ = dim;
    len_ = 1;
    // linspace and prepare possible entries */
    entries_.resize(dim);
    resolutions_.resize(dim);
    step_size_.resize(dim);
    for (int i{}; i < dim; ++i) {
      const int& res = resolutions[i];
      if (res < 2) {
        splinepy::utils::PrintAndThrowError(
            "Resolutions for GridPoints can't be less than "
            "2.");
      }
      len_ *= res;
      resolutions_[i] = res;

      Vector_<double>& entryvec = entries_[i];
      entryvec.resize(res);

      double& step_size = step_size_[i];
      // get current lower and upper bounds
      const double& lower = bounds[i];
      const double& upper = bounds[i + dim];
      // get step size
      step_size = (upper - lower) / static_cast<double>(res - 1);
      // linspace
      for (int j{}; j < res; ++j) {
        entryvec[j] = lower + (step_size * j);
      }
    }
  }

  /// @brief given global id, returns grid point coordinate
  /// @param id in
  /// @param grid_point out
  void IdToGridPoint(const int& id, double* grid_point) const {
    int tmp{id};
    for (int i{0}; i < dim_; ++i) {
      grid_point[i] = entries_[i][tmp % resolutions_[i]];
      tmp /= resolutions_[i];
    }
  }

  /// @brief fills array with full set of grid points
  /// @param grid_point_to_fill
  void Fill(double* grid_point_to_fill) const {
    int i{}, tile{len_}, repeat{1};
    for (const auto& entry : entries_) {
      const int entry_size = entry.size();
      tile /= entry_size;
      int output_id{i};
      for (int j{}; j < tile; ++j) {
        for (const auto& e : entry) {
          for (int k{}; k < repeat; ++k) {
            grid_point_to_fill[output_id] = e;
            output_id += dim_;
          }
        }
      }
      ++i;
      repeat *= entry_size;
    }
  }

  /// @brief Get global IDs
  /// @param grid_resolutions
  /// @param plane_normal_axis
  /// @param plane_id
  static Vector_<int> IdsOnHyperPlane(const int* grid_resolutions,
                                      const int& dim,
                                      const int& plane_normal_axis,
                                      const int& plane_id) {
    // Determine size of return vector
    int n_points_on_plane{1};
    for (int i_pd{}; i_pd < dim; i_pd++) {
      if (i_pd == plane_normal_axis)
        continue;
      n_points_on_plane *= grid_resolutions[i_pd];
    }

    auto local_to_global_bd_id = [&](const int& local_id) {
      int combined_offsets{1}, id{local_id}, global_id{};
      for (int i_pdc{}; i_pdc < dim; ++i_pdc) {
        const int& res = grid_resolutions[i_pdc];
        if (i_pdc == plane_normal_axis) {
          global_id += combined_offsets * plane_id;
        } else {
          // Determine index in coordinate indexing system
          const int rasterized_index = id % res;
          // Use the current id to update the global index
          global_id += rasterized_index * combined_offsets;
          // Update the copy of the index to search for next slice
          id -= rasterized_index;
          id /= res;
        }
        combined_offsets *= res;
      }
      return global_id;
    };

    // Fill return list
    Vector_<int> return_list(n_points_on_plane);
    for (int i{}; i < n_points_on_plane; ++i) {
      return_list[i] = local_to_global_bd_id(i);
    }
    return return_list;
  }

  /// @brief Size
  int Size() const { return len_; }
};

} /* namespace splinepy::utils */
