#include "splinepy/utils/coordinate_pointers.hpp"

#include "splinepy/utils/print.hpp"

namespace splinepy::utils {

int ControlPointPointers::Len() const { return coordinate_begins_.size(); }

int ControlPointPointers::Dim() const {
  assert(dim_ > 0);
  return dim_;
}

void ControlPointPointers::SetRow(const int id, const double* values) {
  if (invalid_) {
    return;
  }

  if (for_rational_) {
    const double& weight = *(weight_pointers_->weights_[id]);
    double* coord = coordinate_begins_[id];
    for (int i{}; i < Dim(); ++i) {
      coord[i] = values[i] * weight;
    }
  } else {
    double* coord = coordinate_begins_[id];
    for (int i{}; i < Dim(); ++i) {
      coord[i] = values[i];
    }
  }
}

void ControlPointPointers::Sync(const double* values) {
  if (invalid_) {
    return;
  }

  const int dim = Dim();

  if (for_rational_) {
    const auto& weights = weight_pointers_->weights_;

    for (int i{}; i < Len(); ++i) {
      // get destinations and sources
      double* current_coord = coordinate_begins_[i];
      const double* current_value = &values[i * dim];
      const double& current_weight = *weights[i];

      for (int j{}; j < dim; ++j) {
        // saves weighted control points
        current_coord[j] = current_value[j] * current_weight;
      }
    }
  } else {
    for (int i{}; i < Len(); ++i) {
      // get destinations and sources
      double* current_coord = coordinate_begins_[i];
      const double* current_value = &values[i * dim];

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

void WeightPointers::SetRow(const int id, double const& value) {
  if (invalid_) {
    return;
  }

  if (auto cpp = control_point_pointers_.lock()) {
    // adjustment factor - new value divided by previous factor;
    double& current_weight = *weights_[id];
    const double adjust_factor = value / current_weight;

    double* current_coordinate = cpp->coordinate_begins_[id];
    for (int i{}; i < cpp->dim_; ++i) {
      current_coordinate[i] *= adjust_factor;
    }

    // save new weight
    current_weight = value;
  } else {
    splinepy::utils::PrintAndThrowError(
        "Missing related control point pointers. Please help us and report "
        "this issue to github.com/tataratat/splinepy, thank you!");
  }
}

void WeightPointers::Sync(const double* values) {
  if (invalid_) {
    return;
  }
  for (int i{}; i < Len(); ++i) {
    SetRow(i, values[i]);
  }
}

} // namespace splinepy::utils
