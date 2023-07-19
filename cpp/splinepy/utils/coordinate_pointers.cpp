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

std::shared_ptr<ControlPointPointers>
ControlPointPointers::SubSetIncomplete(const int* ids, const int n_ids) {
  auto subset = std::make_shared<ControlPointPointers>();
  subset->coordinate_begins_.reserve(n_ids);
  subset->dim_ = dim_;
  subset->for_rational_ = for_rational_;

  for (int i{}; i < n_ids; ++i) {
    subset->coordinate_begins_.push_back(coordinate_begins_[ids[i]]);
  }
  return subset;
}

std::shared_ptr<ControlPointPointers>
ControlPointPointers::SubSet(const int* ids, const int n_ids) {
  auto subset = SubSetIncomplete(ids, n_ids);
  if (for_rational_) {
    subset->weight_pointers_ = weight_pointers_->SubSetIncomplete(ids, n_ids);
    weight_pointers_->control_point_pointers_ = subset;
  }
  return subset;
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

std::shared_ptr<WeightPointers>
WeightPointers::SubSetIncomplete(const int* ids, const int n_ids) {
  auto subset = std::make_shared<WeightPointers>();
  subset->weights_.reserve(n_ids);

  for (int i{}; i < n_ids; ++i) {
    subset->weights_.push_back(weights_[ids[i]]);
  }
  return subset;
}

std::shared_ptr<WeightPointers> WeightPointers::SubSet(const int* ids,
                                                       const int n_ids) {
  auto subset = SubSetIncomplete(ids, n_ids);
  subset->control_point_pointers_ =
      control_point_pointers_->SubSetIncomplete(ids, n_ids);
  subset->control_point_pointers_->weight_pointers_ = subset;
  return subset;
}

} // namespace splinepy::utils
