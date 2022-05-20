#pragma once

#include <cmath>
#include <thread>
#include <utility>


/* thread for all */
template <typename Func, typename IndexT>
void nthread_execution(
    Func& f,
    const IndexT& total,
    const IndexT& nthread
) {
  // if nthread == 1, don't even bother creating thread
  if (nthread == 1) {
    f(0, total);
    return;
  }

  // get chunk size and prepare threads
  const IndexT chunk_size = std::ceil(total / nthread);
  std::vector<std::thread> tpool;
  tpool.reserve(nthread);

  for (int i{0}; i < (nthread - 1); i++) {
    tpool.emplace_back(std::thread{f, i * chunk_size, (i + 1)*chunk_size});
  }
  {
    // last one
    tpool.emplace_back(
        std::thread{f, (nthread - 2) * chunk_size, total}
    );
  }

  for (auto &t : tpool) {
    t.join();
  }
}

/* raster point generator */
template <int dim, typename FillType>
void raster_points(
    
) {
}


template <typename DataT,
          typename IndexT,
          int dim>
class RasterPoints {
public:
  RasterPoints() = default;
  RasterPoints(std::array<std::array<DataT, dim>, 2>& bounds,
               std::array<IndexT, dim>& resolutions) : res_(resolutions){
    // linspace and prepare possible entries */
    len_ = 1;
    for (int i{0}; i < dim; i++) {
      len_ *= resolutions[i];
      std::vector<DataT> entryvec{};
      entryvec.reserve(resolutions[i]);

      step_size_[i] = (bounds[1][i] - bounds[0][i]) / (res_[i] - 1);
      for (int j{0}; j < resolutions[i]; j++) {
        entryvec.emplace_back(bounds[0][i] + step_size_[i] * j);
      }
      entries_[i] = std::move(entryvec);
    }
  }


  const std::array<DataT, dim> operator[](IndexT id) {
    IndexT quot = id, newquot, rem;
    std::array<DataT, dim> out;

    for (int i{0}; i < dim; i++) {
      // compiler should be smart enough
      newquot = (IndexT) quot / res_[i];
      out[i] = entries_[i][quot % res_[i]];
      quot = newquot;
    }

    return std::move(out);
  }

  template<typename ParaCoord>
  void id_to_paracoord(const IndexT& id, ParaCoord& pcoord) {
    using SPC = typename ParaCoord::value_type;

    IndexT quot = id, newquot, rem;
    for (int i{0}; i < dim; i++) {
      newquot = (IndexT) quot / res_[i];
      pcoord[i] = SPC{entries_[i][quot % res_[i]]};
      quot = newquot;
    }
  }

  IndexT size() {
    return len_;
  }

  IndexT len() {
    return len_;
  }



private:
  std::array<DataT, dim> step_size_; /* don't really need it .. */
  std::array<IndexT, dim> res_;
  IndexT count_;
  IndexT len_;
  std::array<std::vector<DataT>, dim> entries_;

};
