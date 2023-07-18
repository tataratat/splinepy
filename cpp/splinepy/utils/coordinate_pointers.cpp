#include "splinepy/utils/coordinate_pointers.hpp"

namespace splinepy::utils {

int ControlPointPointers::Len() const { return coordinate_begins_.size(); }
int ControlPointPointers::Dim() const {
  assert(dim_ > 0);
  return dim_;
}

void ControlPointPointers::SetRow(const int id, const double* values) {
  if (for_rational_) {
    const auto& weight = *(weight_pointers_->weights_[id]);
    auto* coord = coordinate_begins_[id];
    for (int i{}; i < Dim(); ++i) {
      coord[i] = values[i] * weight;
    }
  } else {
    auto* coord = coordinate_begins_[id];
    for (int i{}; i < Dim(); ++i) {
      coord[i] = values[i];
    }
  }
}

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
      const auto* current_value = &values[current_id * dim];
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
      const auto* current_value = &values[current_id * dim];

      for (int j{}; j < dim; ++j) {
        // saves weighted control points
        current_coord[j] = current_value[j];
      }
    }
  }
}

void ControlPointPointers::Sync(const double* values) {
  const auto dim = Dim();

  if (for_rational_) {
    const auto& weights = weight_pointers_->weights_;

    for (int i{}; i < Len(); ++i) {
      // get destinations and sources
      auto* current_coord = coordinate_begins_[i];
      const auto* current_value = &values[i * dim];
      const auto& current_weight = *weights[i];

      for (int j{}; j < dim; ++j) {
        // saves weighted control points
        current_coord[j] = current_value[j] * current_weight;
      }
    }
  } else {
    for (int i{}; i < Len(); ++i) {
      // get destinations and sources
      auto* current_coord = coordinate_begins_[i];
      const auto* current_value = &values[i * dim];

      for (int j{}; j < dim; ++j) {
        // saves weighted control points
        current_coord[j] = current_value[j];
      }
    }
  }
}

int WeightPointers::Len() const { return weights_.size(); }
int WeightPointers::Dim() const {
  assert(dim_ > 0);
  return dim_;
}

void WeightPointers::SetRow(const int id, const double value) {
  // adjustment factor - new value divided by previous factor;
  auto& current_weight = *weights_[id];
  const double adjust_factor = value / current_weight;

  double* current_coordinate = control_point_pointers_->coordinate_begins_[id];
  for (int i{}; i < control_point_pointers_->dim_; ++i) {
    current_coordinate[i] *= adjust_factor;
  }

  // save new weight
  current_weight = value;
}

void WeightPointers::SetRows(const int* ids,
                             const int n_rows,
                             const double* values) {
  for (int i{}; i < n_rows; ++i) {
    SetRow(ids[i], values[i]);
  }
}

void WeightPointers::Sync(const double* values) {
  for (int i{}; i < Len(); ++i) {
    SetRow(i, values[i]);
  }
}

} // namespace splinepy::utils
