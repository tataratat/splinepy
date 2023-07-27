#pragma once

#include <cmath>
#include <numeric>
#include <utility>

#include <splinepy/utils/print.hpp>

namespace splinepy::utils {

/// @brief Grid points class
/// @tparam DataT
/// @tparam IndexType
/// @tparam dim
template<typename DataT, typename IndexType, int dim>
class GridPoints {
public:
  /// @brief Step size, could be useful for setting searchbounds aggressively
  std::array<DataT, dim> step_size_;

  GridPoints() = default;
  /// @brief Constructor
  /// @param bounds
  /// @param resolutions
  GridPoints(const std::array<std::array<DataT, dim>, 2>& bounds,
             const std::array<IndexType, dim>& resolutions)
      : res_(resolutions),
        bounds_(bounds) {
    // linspace and prepare possible entries */
    len_ = 1;
    for (int i{0}; i < dim; i++) {
      len_ *= resolutions[i];
      std::vector<DataT>& entryvec = entries_[i];
      entryvec.reserve(resolutions[i]);

      step_size_[i] = (bounds[1][i] - bounds[0][i]) / (res_[i] - 1);
      for (int j{0}; j < resolutions[i]; j++) {
        entryvec.emplace_back(bounds[0][i] + step_size_[i] * j);
      }
    }
    len_times_dim_ = len_ * dim;
  }

  /// @brief Get grid point from index
  const std::array<DataT, dim> operator[](const IndexType id) {
    return IndexToGridPoint<std::array<DataT, dim>>(id);
  }

  /// @brief Get grid point from index
  /// @param id
  /// @return Grid point
  template<typename GridPointType>
  GridPointType IndexToGridPoint(const IndexType& id) const {

    using ValueType = typename GridPointType::value_type;

    GridPointType gpoint{};

    IndexType quot{id};
    for (int i{0}; i < dim; i++) {
      gpoint[i] = ValueType{entries_[i][quot % res_[i]]};
      quot /= res_[i];
    }

    return gpoint;
  }

  /// subroutine variation
  template<typename GridPointType>
  void IndexToGridPoint(const IndexType& id, GridPointType& gpoint) const {

    using ValueType = typename GridPointType::value_type;

    IndexType quot{id};
    for (int i{0}; i < dim; i++) {
      gpoint[i] = ValueType{entries_[i][quot % res_[i]]};
      quot /= res_[i];
    }
  }

  /// Saved at the first call then returns
  const std::vector<DataT>& GetAllGridPoints() {
    // early exit if we have them
    if (saved_grid_points_.size() == len_times_dim_) {
      return saved_grid_points_;
    }
    // get all.
    saved_grid_points_.clear();
    saved_grid_points_.reserve(len_times_dim_);
    for (IndexType i{}; i < len_; ++i) {
      // third time is the charm
      IndexType quot{i};
      for (IndexType j{}; j < dim; ++j) {
        const auto& r = res_[j];
        saved_grid_points_.emplace_back(entries_[j][quot % r]);
        quot /= r;
      }
    }

    return saved_grid_points_;
  }

  /// @brief Get global IDs
  /// @param grid_resolutions
  /// @param plane_normal_axis
  /// @param plane_id
  static std::vector<IndexType>
  GridPointIdsOnHyperPlane(const std::array<IndexType, dim>& grid_resolutions,
                           const IndexType& plane_normal_axis,
                           const IndexType& plane_id) {
    // Determine size of return vector
    IndexType n_points_on_boundary{1};
    for (std::size_t i_pd{}; i_pd < dim; i_pd++) {
      if (i_pd == plane_normal_axis)
        continue;
      n_points_on_boundary *= grid_resolutions[i_pd];
    }

    auto local_to_global_bd_id = [&](const IndexType& local_id) {
      IndexType combined_offsets{1}, id{local_id}, global_id{};
      for (std::size_t i_pdc{}; i_pdc < dim; ++i_pdc) {
        const auto& res = grid_resolutions[i_pdc];
        if (i_pdc == plane_normal_axis) {
          global_id += combined_offsets * plane_id;
        } else {
          // Determine index in coordinate indexing system
          const IndexType rasterized_index = id % res;
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
    std::vector<IndexType> return_list{};
    return_list.reserve(n_points_on_boundary);
    for (std::size_t i{}; i < n_points_on_boundary; ++i) {
      return_list.push_back(local_to_global_bd_id(i));
    }
    return return_list;
  }

  /// @brief Get IDs of grid points on the boundary
  static std::vector<IndexType>
  GridPointIdsOnBoundary(const std::array<IndexType, dim>& grid_resolutions,
                         const IndexType& plane_normal_axis,
                         const IndexType& extrema) {
    const IndexType plane_id =
        (extrema > 0) ? grid_resolutions[plane_normal_axis] - 1 : 0;
    return GridPointIdsOnHyperPlane(grid_resolutions,
                                    plane_normal_axis,
                                    plane_id);
  }

  /// @brief Size
  IndexType Size() const { return len_; }

  /// @brief Length
  IndexType Len() const { return len_; }

private:
  std::array<IndexType, dim> res_;
  std::array<std::array<DataT, dim>, 2> bounds_;
  IndexType len_;
  IndexType len_times_dim_;
  std::array<std::vector<DataT>, dim> entries_;
  std::vector<DataT> saved_grid_points_;
};

/// @brief CStyleArrayPointer based dynamic grid point sampler
/// Values are kept as int since they come as int from python
class CStyleArrayPointerGridPoints {
public:
  CStyleArrayPointerGridPoints() = default;
  CStyleArrayPointerGridPoints(const int dim,
                               const double* bounds,
                               const int* resolutions) {
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
    entries_.clear();
    entries_.assign(dim, std::vector<double>{});
    resolutions_.clear();
    resolutions_.reserve(dim);
    for (int i{}; i < dim; ++i) {
      const int& res = resolutions[i];
      if (res < 2) {
        splinepy::utils::PrintAndThrowError(
            "Resolutions for CStyleArrayPointerGridPoints can't be less than "
            "2.");
      }
      len_ *= res;
      resolutions_.push_back(res);

      std::vector<double>& entryvec = entries_[i];
      entryvec.reserve(res);

      const double step_size = (bounds[dim + i] - bounds[i]) / (res - 1);
      for (int j{}; j < resolutions[i]; ++j) {
        entryvec.emplace_back(bounds[i] + step_size * j);
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

  /// @brief Size
  int Size() const { return len_; }

protected:
  /// @brief Vector of bounds
  std::vector<double> bounds_;
  /// @brief Vector of resolutions
  std::vector<int> resolutions_;
  /// @brief Vector of vectors of entries
  std::vector<std::vector<double>> entries_;
  /// @brief Length
  int len_;
  /// @brief Dimension
  int dim_;
};

} /* namespace splinepy::utils */
