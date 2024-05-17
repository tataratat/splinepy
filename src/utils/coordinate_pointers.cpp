/*
MIT License

Copyright (c) 2021 Jaewook Lee

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "splinepy/utils/coordinate_pointers.hpp"

namespace splinepy::utils {

void ControlPointPointers::SetRow(const int id_, const double* values) {
  if (invalid_) {
    return;
  }

  // allow negative id
  const int len = Len();
  const int id = WrapId(id_, len);

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

void WeightPointers::SetRow(const int id_, double const& value) {
  if (invalid_) {
    return;
  }

  // allow negative id
  const int len = Len();
  const int id = WrapId(id_, len);

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
