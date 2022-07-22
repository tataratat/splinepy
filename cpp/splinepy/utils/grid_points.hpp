#pragma once

#include <cmath>
#include <utility>


namespace splinepy::utils {

template <typename DataT,
          typename IndexT,
          int dim>
class GridPoints {
public:

  // could be useful for setting searchbounds aggresively
  std::array<DataT, dim> step_size_;

  GridPoints() = default;
  GridPoints(const std::array<std::array<DataT, dim>, 2>& bounds,
             const std::array<IndexT, dim>& resolutions) : res_(resolutions){
    // linspace and prepare possible entries */
    len_ = 1;
    for (int i{0}; i < dim; i++) {
      len_ *= resolutions[i];
      std::vector<DataT>& entryvec = entries[i];
      entryvec.reserve(resolutions[i]);

      step_size_[i] = (bounds[1][i] - bounds[0][i]) / (res_[i] - 1);
      for (int j{0}; j < resolutions[i]; j++) {
        entryvec.emplace_back(bounds[0][i] + step_size_[i] * j);
      }
    }
  }


  const std::array<DataT, dim> operator[](const IndexT id) {
    return IndexToParametricCoordinate<std::array<DataT, dim>>(id);
  }

  template<typename ParaCoord>
  ParaCoord IndexToParametricCoordinate(const IndexT& id) const {

    using ValueType = typename ParaCoord::value_type;

    ParaCoord pcoord{};

    IndexT quot{id};
    for (int i{0}; i < dim; i++) {
      pcoord[i] = ValueType{entries_[i][quot % res_[i]]};
      quot /= res_[i];
    }

    return pcoord;
  }

  /// subroutine variation
  template<typename ParaCoord>
  void IndexToParametricCoordinate(const IndexT& id,
                                   ParaCoord& pcoord) const {

    using ValueType = typename ParaCoord::value_type;

    IndexT quot{id};
    for (int i{0}; i < dim; i++) {
      pcoord[i] = ValueType{entries_[i][quot % res_[i]]};
      quot /= res_[i];
    }

  }


  IndexT Size() const {
    return len_;
  }

  IndexT Len() const {
    return len_;
  }



private:
  std::array<IndexT, dim> res_;
  IndexT count_;
  IndexT len_;
  std::array<std::vector<DataT>, dim> entries_;

};

} /* namespace splinepy::utils */
