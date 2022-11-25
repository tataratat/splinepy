#pragma once

#include <cmath>
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
      : res_(resolutions) {
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

  IndexT Size() const { return len_; }

  IndexT Len() const { return len_; }

private:
  std::array<IndexT, dim> res_;
  IndexT count_;
  IndexT len_;
  std::array<std::vector<DataT>, dim> entries_;
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
