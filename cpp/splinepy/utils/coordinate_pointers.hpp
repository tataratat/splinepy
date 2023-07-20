#pragma once

#include <cassert>
#include <memory>
#include <vector>

namespace splinepy::utils {

struct WeightPointers;

struct ControlPointPointers {
  std::vector<double*> coordinate_begins_;
  int dim_{-1};
  bool for_rational_{false};
  std::shared_ptr<WeightPointers> weight_pointers_ = nullptr;

  int Len() const;
  int Dim() const;
  void SetRow(const int id, const double* values);
  template<bool same_sized_value>
  void SetRows(const int* ids, const int n_rows, const double* values);
  void Sync(const double* values);

  // internal use
  std::shared_ptr<ControlPointPointers> SubSetIncomplete(const int* ids,
                                                         const int n_ids);
  std::shared_ptr<ControlPointPointers> SubSet(const int* ids, const int n_ids);
};

struct WeightPointers {
  std::vector<double*> weights_;
  std::shared_ptr<ControlPointPointers> control_point_pointers_ = nullptr;
  static const int dim_{1};

  int Len() const;
  int Dim() const;
  void SetRow(const int id, double const& value);
  template<bool same_sized_value>
  void SetRows(const int* ids, const int n_rows, const double* values);
  void Sync(const double* values);
  // internal use
  std::shared_ptr<WeightPointers> SubSetIncomplete(const int* ids,
                                                   const int n_ids);
  std::shared_ptr<WeightPointers> SubSet(const int* ids, const int n_ids);
};

template<bool same_sized_values>
void ControlPointPointers::SetRows(const int* ids,
                                   const int n_rows,
                                   const double* values) {
  const auto dim = Dim();

  if (for_rational_) {
    const auto& weights = weight_pointers_->weights_;

    for (int i{}; i < n_rows; ++i) {
      // get id
      const auto& current_id = ids[i];
      // get destinations and sources
      auto* current_coord = coordinate_begins_[current_id];
      const double* current_value;
      if constexpr (same_sized_values) {
        current_value = &values[current_id * dim];
      } else {
        current_value = &values[i * dim];
      }
      const auto& current_weight = *weights[current_id];

      for (int j{}; j < dim; ++j) {
        // saves weighted control points
        current_coord[j] = current_value[j] * current_weight;
      }
    }
  } else {
    for (int i{}; i < n_rows; ++i) {
      // get id
      const auto& current_id = ids[i];
      // get destinations and sources
      auto* current_coord = coordinate_begins_[current_id];
      const double* current_value;
      if constexpr (same_sized_values) {
        current_value = &values[current_id * dim];
      } else {
        current_value = &values[i * dim];
      }

      for (int j{}; j < dim; ++j) {
        // saves weighted control points
        current_coord[j] = current_value[j];
      }
    }
  }
}

template<bool same_sized_values>
void WeightPointers::SetRows(const int* ids,
                             const int n_rows,
                             const double* values) {
  for (int i{}; i < n_rows; ++i) {
    if constexpr (same_sized_values) {
      SetRow(ids[i], values[ids[i]]);
    } else {
      SetRow(ids[i], values[i]);
    }
  }
}

} // namespace splinepy::utils
