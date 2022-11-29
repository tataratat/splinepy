#pragma once

#include <cmath>
#include <numeric>
#include <utility>

#include <splinepy/utils/print.hpp>

namespace splinepy::utils {

template<typename DataT, typename IndexT, int dim>
class GridPoints {
public:
  // could be useful for setting searchbounds aggresively
  std::array<DataT, dim> step_size_;

  GridPoints() = default;
  GridPoints(const std::array<std::array<DataT, dim>, 2>& bounds,
             const std::array<IndexT, dim>& resolutions)
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

  const std::array<DataT, dim> operator[](const IndexT id) {
    return IndexToGridPoint<std::array<DataT, dim>>(id);
  }

  template<typename GridPointType>
  GridPointType IndexToGridPoint(const IndexT& id) const {

    using ValueType = typename GridPointType::value_type;

    GridPointType gpoint{};

    IndexT quot{id};
    for (int i{0}; i < dim; i++) {
      gpoint[i] = ValueType{entries_[i][quot % res_[i]]};
      quot /= res_[i];
    }

    return gpoint;
  }

  /// subroutine variation
  template<typename GridPointType>
  void IndexToGridPoint(const IndexT& id, GridPointType& gpoint) const {

    using ValueType = typename GridPointType::value_type;

    IndexT quot{id};
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
    for (IndexT i{}; i < len_; ++i) {
      // third time is the charm
      IndexT quot{i};
      for (IndexT j{}; j < dim; ++j) {
        const auto& r = res_[j];
        saved_grid_points_.emplace_back(entries_[j][quot % r]);
        quot /= r;
      }
    }

    return saved_grid_points_;
  }

  static std::vector<IndexT>
  GridPointIdsOnHyperPlane(const std::array<IndexT, dim>& grid_resolutions,
                           const IndexT& plane_normal_axis,
                           const IndexT& plane_id) {
    // Determine size of return vector
    IndexT n_points_on_boundary{1};
    for (std::size_t i_pd{}; i_pd < dim; i_pd++) {
      if (i_pd == plane_normal_axis)
        continue;
      n_points_on_boundary *= grid_resolutions[i_pd];
    }

    auto local_to_global_bd_id = [&](const IndexT& local_id) {
      IndexT offset{1}, id{local_id}, global_id{};
      for (std::size_t i_pdc{}; i_pdc < dim; ++i_pdc) {
        const auto& res = grid_resolutions[i_pdc];
        if (i_pdc == plane_normal_axis) {
          global_id += offset * plane_id;
        } else {
          const IndexT i = id % res;
          global_id += id * offset;
          id -= i;
          id /= res;
        }
        offset *= res;
      }
      return global_id;
    };

    // Fill return list
    std::vector<IndexT> return_list{};
    return_list.reserve(n_points_on_boundary);
    for (std::size_t i{}; i < n_points_on_boundary; ++i) {
      return_list.push_back(local_to_global_bd_id(i));
    }
    return return_list;
  }

  /// Extract grid point ids that lies on specified hyper plane
  std::vector<IndexT> GridPointIdsOnHyperPlane(const IndexT& plane_normal_axis,
                                               const IndexT& plane_id) {
    IndexT n_plane_points{1};
    for (IndexT i{}; i < dim; ++i) {
      if (plane_normal_axis == i) {
        continue;
      }
      n_plane_points *= res_[i];
    }

    // prepare return
    std::vector<IndexT> hits;
    hits.reserve(n_plane_points);

    // get search value
    const DataT& search_value = entries_[plane_normal_axis][plane_id];

    const auto& all_grid_points = GetAllGridPoints();

    // loop and see if value matches, if so, append.
    for (IndexT i{}; i < len_; ++i) {
      if constexpr (std::numeric_limits<DataT>::is_integer) {
        if (all_grid_points[i * dim + plane_normal_axis] == search_value) {
          hits.emplace_back(i);
        }
      } else {
        if (std::abs(all_grid_points[i * dim + plane_normal_axis]
                     - search_value)
            < 1e-12) {
          hits.emplace_back(i);
        }
      }
    }

    // size check
    if (hits.size() != n_plane_points) {
      splinepy::utils::PrintAndThrowError(
          "Sorry, something went wrong during hyper plane id extraction. "
          "Expecting",
          n_plane_points,
          "but only found",
          hits.size());
    }

    // all good
    return hits;
  }

  /// Extract boundary Ids
  std::vector<IndexT> GridPointIdsOnBoundary(const IndexT& plane_normal_axis,
                                             const IndexT& extrema) {
    // findout which bounds - this matches proximity clip description.
    // in case IndexT is unsigned, we check if extrema is greater than zero.
    // if yes, upper, else, lower
    const IndexT plane_id =
        (extrema > 0) ? entries_[plane_normal_axis].size() - 1 : 0;
    return GridPointIdsOnHyperPlane(plane_normal_axis, plane_id);
  }

  IndexT Size() const { return len_; }

  IndexT Len() const { return len_; }

private:
  std::array<IndexT, dim> res_;
  std::array<std::array<DataT, dim>, 2> bounds_;
  IndexT len_;
  IndexT len_times_dim_;
  std::array<std::vector<DataT>, dim> entries_;
  std::vector<DataT> saved_grid_points_;
};

/// CStyleArrayPointer based dynamic grid point sampler
/// Values are kept as int since they come as int from python
class CStyleArrayPointerGridPoints {
public:
  CStyleArrayPointerGridPoints() = default;
  CStyleArrayPointerGridPoints(const int dim,
                               const double* bounds,
                               const int* resolutions) {
    SetUp(dim, bounds, resolutions);
  }

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

  void IdToGridPoint(const int& id, double* grid_point) const {
    int tmp{id};
    for (int i{0}; i < dim_; ++i) {
      grid_point[i] = entries_[i][tmp % resolutions_[i]];
      tmp /= resolutions_[i];
    }
  }

  int Size() const { return len_; }

protected:
  std::vector<double> bounds_;
  std::vector<int> resolutions_;
  std::vector<std::vector<double>> entries_;
  int len_;
  int dim_;
};

} /* namespace splinepy::utils */
