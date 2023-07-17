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
  void SetRows(const int* ids, const int n_rows, const double* values);
  void Sync(const double* values);
};

struct WeightPointers {
  std::vector<double*> weights_;
  std::shared_ptr<ControlPointPointers> control_point_pointers_ = nullptr;
  static const int dim_{1};

  int Len() const;
  int Dim() const;
  void SetRow(const int id, const double value);
  void SetRows(const int* ids, const int n_rows, const double* values);
  void Sync(const double* values);
};

} // namespace splinepy::utils
