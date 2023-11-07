#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <type_traits>

#include "splinepy/utils/print.hpp"

namespace splinepy::utils {

/// @brief fully dynamic array
/// @tparam DataType
/// @tparam dim
template<typename DataType, int dim = 1, typename IndexType = int>
class Array {
  static_assert(dim > 0, "dim needs to be positive value bigger than zero.");
  static_assert(std::is_integral_v<IndexType>,
                "IndexType should be an integral type");

public:
  using ShapeType_ = std::array<IndexType, dim>;
  using StridesType_ = std::array<IndexType, dim - 1>;
  static constexpr const IndexType kDim = static_cast<IndexType>(dim);

  using DataType_ = DataType;
  using IndexType_ = IndexType;

  // std container like types
  using value_type = DataType;
  using size_type = IndexType;

protected:
  bool own_data_{false};
  /// @brief beginning of the array pointer that this object manages
  DataType* data_{nullptr};

  /// @brief size of this Array
  IndexType size_{};

  /// @brief strides in case this is a higher dim. last entry should be the same
  /// as size_
  StridesType_ strides_;

  /// @brief shape of the array
  ShapeType_ shape_;

  /// @brief helper to find id offset
  /// @tparam ...Ts
  /// @param index
  /// @param counter
  /// @param ...id
  template<typename... Ts>
  constexpr void ComputeRemainingOffset(IndexType& index,
                                        IndexType& counter,
                                        const IndexType& id0,
                                        const Ts&... id) const {
    // last element should be just added
    if constexpr (sizeof...(Ts) == 0) {
      index += id0;
      // intentionally skipping
      // ++counter;
      return;
    }

    // add strides * id
    index += strides_[counter++] * id0;

    // do the rest
    if constexpr (sizeof...(Ts) > 0) {
      ComputeRemainingOffset(index, counter, id...);
    }
  }

public:
  constexpr DataType* data() {
    assert(data_);
    return data_;
  }

  constexpr const DataType* data() const {
    assert(data_);
    return data_;
  }

  constexpr DataType* begin() const {
    assert(data_);
    return data_;
  }

  constexpr DataType* end() const {
    assert(data_);
    assert(size_ > 0);
    return data_ + size_;
  }

  constexpr IndexType size() const { return size_; }

  constexpr void DestroyData() {
    if (own_data_ && data_) {
      delete[] data_;
    }

    data_ = nullptr;
    own_data_ = false;
  }

  constexpr void SetData(DataType* data_pointer) {
    DestroyData();
    data_ = data_pointer;
  }

  constexpr DataType* GetData() {
    assert(data_);
    return data_;
  }

  constexpr const DataType* GetData() const {
    assert(data_);
    return data_;
  }

  /// @brief reallocates space. After this call, own_data_ should be true
  /// @param size
  constexpr void Reallocate(const IndexType& size) {
    // destroy and reallocate space
    DestroyData();
    data_ = new DataType[size];
    own_data_ = true;

    // set size - don't forget to set shape in case this is multi-dim array
    size_ = size;

    // set shape if this is a 1d array
    if constexpr (dim == 1) {
      shape_[0] = size;
    }
  }

  template<typename... Ts>
  constexpr void SetShape(const IndexType& shape0, const Ts&... shape) {
    static_assert(sizeof...(Ts) == static_cast<std::size_t>(dim - 1),
                  "shape should match the dim");

    // early exit in case of dim 1
    if constexpr (sizeof...(Ts) == 0) {
      size_ = shape0;
      shape_[0] = shape0;
      return;
    }

    // (re)initialize shape_
    shape_ = {shape0, shape...};

    // (re)initialize  strides_
    strides_ = {shape...};

    // reverse iter and accumulate
    IndexType sub_total{1};
    typename StridesType_::reverse_iterator r_stride_iter = strides_.rbegin();
    for (; r_stride_iter != strides_.rend(); ++r_stride_iter) {
      // accumulate
      sub_total *= *r_stride_iter;
      // assign
      *r_stride_iter = sub_total;
    }

    // at this point, sub total has accumulated shape values except shape0
    // set size
    size_ = sub_total * shape0;
  }

  /// @brief default. use SetData() and SetShape<false>
  Array() = default;

  /// @brief basic array ctor
  /// @tparam ...Ts
  /// @param ...shape
  template<typename... Ts>
  Array(const Ts&... shape) {
    SetShape(shape...);
    Reallocate(size_);
  }

  /// @brief data wrapping array
  /// @tparam ...Ts
  /// @param data_pointer
  /// @param ...shape
  template<typename... Ts>
  Array(DataType* data_pointer, const Ts&... shape) {
    SetData(data_pointer);
    SetShape(shape...);
  }

  /// @brief copy ctor
  /// @param other
  Array(const Array& other) {
    // memory alloc
    Reallocate(other.size_);
    // copy data
    std::copy_n(other.data(), other.size_, data_);
    // copy shape
    shape_ = other.shape_;
    // copy strides
    strides_ = other.strides_;
  }

  /// @brief move ctor
  /// @param other
  Array(Array&& other) {
    own_data_ = other.own_data_;
    data_ = std::move(other.data_);
    size_ = std::move(other.size_);
    strides_ = std::move(other.strides_);
    shape_ = std::move(other.shape_);
    other.own_data_ = false;
  }

  /// @brief copy assignment - need to make sure that you have enough space
  /// @param rhs
  /// @return
  constexpr Array& operator=(const Array& rhs) {
    // size check is crucial, so runtime check
    if (size_ != rhs.size_) {
      splinepy::utils::PrintAndThrowError("Array::operator=(const Array&)",
                                          "size mismatch between rhs");
    }
    std::copy_n(rhs.data(), rhs.size(), data_);

    // copy shape and stride info
    shape_ = rhs.Shape();
    strides_ = rhs.Strides();

    return *this;
  }

  /// @brief move assignment. currently same as ctor.
  /// @param rhs
  /// @return
  constexpr Array& operator=(Array&& rhs) {
    DestroyData();

    own_data_ = rhs.own_data_;
    data_ = std::move(rhs.data_);
    size_ = std::move(rhs.size_);
    strides_ = std::move(rhs.strides_);
    shape_ = std::move(rhs.shape_);
    rhs.own_data_ = false;
    return *this;
  }

  /// @brief
  ~Array() { DestroyData(); }

  /// @brief Returns const shape object
  /// @return
  constexpr const ShapeType_& Shape() const { return shape_; }

  /// @brief Returns const strides object
  /// @return
  constexpr const StridesType_& Strides() const { return strides_; }

  /// @brief flat indexed array access
  /// @param index
  /// @return
  constexpr DataType& operator[](const IndexType& index) {
    assert(index > -1 && index < size_);

    return data_[index];
  }

  /// @brief flat indexed array access
  /// @param index
  /// @return
  constexpr const DataType& operator[](const IndexType& index) const {
    assert(index > -1 && index < size_);

    return data_[index];
  }

  /// @brief multi-indexed array access
  /// @tparam ...Ts
  /// @param id0
  /// @param ...id
  /// @return
  template<typename... Ts>
  constexpr DataType& operator()(const IndexType& id0, const Ts&... id) {
    static_assert(sizeof...(Ts) == static_cast<std::size_t>(dim - 1),
                  "number of index should match dim");

    if constexpr (dim == 1) {
      return operator[](id0);
    }

    IndexType final_id{id0 * strides_[0]}, i{1};
    if constexpr (sizeof...(Ts) > 0) {
      ComputeRemainingOffset(final_id, i, id...);
    }

    assert(final_id > -1 && final_id < size_);

    return data_[final_id];
  }

  /// @brief multi-indexed array access
  /// @tparam ...Ts
  /// @param id0
  /// @param ...id
  /// @return
  template<typename... Ts>
  constexpr const DataType& operator()(const IndexType& id0,
                                       const Ts&... id) const {
    static_assert(sizeof...(Ts) == static_cast<std::size_t>(dim - 1),
                  "number of index should match dim");

    if constexpr (dim == 1) {
      return operator[](id0);
    }

    IndexType final_id{id0 * strides_[0]}, i{1};
    if constexpr (sizeof...(Ts) > 0) {

      ComputeRemainingOffset(final_id, i, id...);
    }

    assert(final_id > -1 && final_id < size_);

    return data_[final_id];
  }

  /// @brief this[:] = v
  /// @param v
  constexpr void Fill(const DataType& v) {
    assert(data_);

    std::fill_n(data_, size_, v);
  }

  /// @brief this[i] += a[i]
  /// @tparam Iterable
  /// @param a
  template<typename Iterable>
  constexpr void Add(const Iterable& a) {
    assert(a.size() == size_);
    assert(data_);

    for (IndexType i{}; i < size_; ++i) {
      data_[i] += a[i];
    }
  }

  /// @brief this[i] -= a[i]
  /// @tparam Iterable
  /// @param a
  template<typename Iterable>
  constexpr void Subtract(const Iterable& a) {
    assert(a.size() == size_);
    assert(data_);

    for (IndexType i{}; i < size_; ++i) {
      data_[i] -= a[i];
    }
  }

  /// @brief c = (this) dot (a).
  /// @tparam Iterable
  /// @param a
  /// @return
  template<typename Iterable>
  constexpr DataType InnerProduct(const Iterable& a) const {
    static_assert(dim == 1,
                  "inner product is only applicable for dim=1 Array.");
    assert(a.size() == size_);
    assert(data_);

    DataType dot{};
    for (IndexType i{}; i < size_; ++i) {
      dot += a[i] * data_[i];
    }

    return dot;
  }

  constexpr void AAt(Array& aa_t) const {
    static_assert(dim == 2, "AAt is only applicable for dim=2 array.");

    const auto& height = shape_[0];
    const auto& width = shape_[1];

    for (int i{}; i < height; ++i) {
      // compute upper triangle
      for (int j{i}; j < height; ++j) {

        // element to fill
        DataType& ij = aa_t(i, j);
        // init
        ij = 0.0;
        // get beginning of i and j row
        const DataType* i_ptr = &operator()(i, 0);
        const DataType* j_ptr = &operator()(j, 0);

        // inner product
        for (int k{}; k < width; ++k) {
          ij += i_ptr[k] * j_ptr[k];
        }

        // fill symmetric part - this will write diagonal twice. it's okay.
        aa_t(j, i) = ij;
      }
    }
  }

  /// @brief given upper and lower bounds, clips data values inplace.
  /// @tparam Iterable
  /// @tparam IntegerType
  /// @param lower_bound
  /// @param upper_bound
  /// @param clipped
  template<typename Iterable, typename IntegerType>
  constexpr void Clip(const Iterable& lower_bound,
                      const Iterable& upper_bound,
                      Array<IntegerType, dim, IndexType>& clipped) {
    static_assert(std::is_integral_v<IntegerType>,
                  "clipped array should be an integral type");
    static_assert(std::is_signed_v<IntegerType>,
                  "clipped array should be a signed type");

    assert(lower_bound.size() == size_);
    assert(upper_bound.size() == size_);

    for (IndexType i{}; i < size_; ++i) {

      // deref values to write
      DataType& data = data_[i];
      IntegerType& clip = clipped[i];

      if (const DataType& ub = upper_bound[i]; data > ub) {
        // upper bound check
        clip = 1;
        data = ub;
      } else if (const DataType& lb = lower_bound[i]; data < lb) {
        // lower bound check
        clip = -1;
        data = lb;
      } else {
        // within the bounds
        clip = 0;
      }
    }
  }

  /// @brief L2 norm
  /// @return
  constexpr DataType NormL2() {
    DataType norm{};
    for (IndexType i{}; i < size_; ++i) {
      const DataType& data_i = data_[i];
      norm += data_i * data_i;
    }
    return std::sqrt(norm);
  }

  /// @brief returns number of non (exact) zero elements.
  /// @return
  constexpr IndexType NonZeros() const {
    IndexType n{};
    for (IndexType i{}; i < size_; ++i) {
      if (n != static_cast<DataType>(0)) {
        ++n;
      }
    }
    return n;
  }
};

/// @brief c[i] = a[i] - b[i]
/// @tparam IterableA
/// @tparam IterableB
/// @tparam IterableC
/// @param a
/// @param b
/// @param c
template<typename IterableA, typename IterableB, typename IterableC>
inline void Subtract(const IterableA& a, const IterableB& b, IterableC& c) {
  assert(a.size() == b.size());
  assert(b.size() == c.size());

  const auto len = c.size();

  for (std::remove_const_t<decltype(len)> i{}; i < len; ++i) {
    c[i] = a[i] - b[i];
  }
}

template<typename DataType, typename IndexType>
class Matrix {
  using DataType_ = DataType;
  using IndexType_ = IndexType;
  using SignedIndexType = std::make_signed_t<IndexType>;
  using DataArray_ = Array<DataType, 2, IndexType>;
  using RowOrderArray_ = Array<IndexType, 1, IndexType>;

protected:
  /// @brief 2d base array of this matrix
  DataArray_& array_;

  /// @brief row permutation.
  RowOrderArray_ row_order_;

public:
  Matrix(DataArray_& array) : array_(array) {
    const IndexType& height = array.Shape()[0];
    // alloc n_rows
    row_order_.Reallocate(height);
    // arange n_rows
    ResetRowOrder();
  }

  constexpr void ResetRowOrder() {
    std::iota(row_order_.begin(), row_order_.end(), 0);
  }

  constexpr void SwapRow(const IndexType& i, const IndexType j) {
    std::swap(row_order_[i], row_order_[j]);
  }

  constexpr const RowOrderArray_& RowOrder() const { return row_order_; }

  template<typename IterableB, typename IterableX>
  constexpr void Solve(IterableB& b, IterableX& x) {

    // get length of the array/matrix
    const IndexType& len = array_.Shape()[0];

    // array_ should be a square matrix
    assert(len == array_.Shape()[1]);
    // len should match the size of b and x
    assert(len == b.size());
    assert(len == x.size());

    // partial pivoting and forward reduction
    // since we reorder indices only, we append *_r to reordered ids
    for (IndexType i{0}; i < len; i++) {
      // sneak in x initialization here.
      // it doesn't matter if this is i or i_r
      x[i] = 0.;

      // ignore clipped entries
      // if (skipmask[i] != 0)
      //  continue;

      /* swap */
      // go through the rows and compare magnitude and find max row at col i.
      // we use reordered index here to retrieve rows, but keep the
      // true index for everything else
      const IndexType& i_r = row_order_[i];
      IndexType max_row{i};
      DataType current_max{std::abs(array_(i_r, i))};
      // loop with true row index
      for (IndexType j{i + 1}; j < len; ++j) {
        // but look up reordered row's value
        DataType maybe_max{std::abs(array_(row_order_[j], i))};

        if (maybe_max > current_max) {
          current_max = std::move(maybe_max);
          max_row = j;
        }
      }
      // swap if needed. max entry row goes to i-th row.
      if (max_row != i) {
        SwapRow(max_row, i);
        // std::swap(skipmask[maxrow], skipmask[i]);
      }
      /* END swap */

      // forward reduction
      const DataType a_ii = array_(i_r, i);
      for (IndexType j{i + 1}; j < len; ++j) {
        // retrieved reordered j-th row
        const IndexType& j_r = row_order_[j];
        // reduction factor for this row
        const DataType reduction_factor = array_(j_r, i) / a_ii;

        // matrix reduction
        // get the beginnings of i and j rows
        DataType* a_i = &array_(i_r, 0);
        DataType* a_j = &array_(j_r, 0);
        // start with setting the i th column zero.
        a_j[i] = 0.;
        // this is column entry, so true indices
        // reduce the rest of the column
        for (IndexType k{i + 1}; k < len; ++k) {
          a_j[k] -= a_i[k] * reduction_factor;
        }

        // row-vec reduction
        b[j_r] -= b[i_r] * reduction_factor;
      }
      /* END forward reduction */
    }

    // back substitution
    // Loop backwards
    for (SignedIndexType i{len - 1}; i > static_cast<SignedIndexType>(-1);
         --i) {
      // ignore clipped entries
      // if (skipmask[i] != 0)
      // continue;

      const IndexType& i_r = row_order_[i];
      DataType* a_i = &array_(i_r, 0);
      DataType sum{};
      for (IndexType j{i + 1}; j < len; ++j) {
        sum += a_i[j] * x[row_order_[j]];
      }
      x[i_r] = (b[i_r] - sum) / a_i[i];
    }
  }
};

/// @brief Gauss elimination with partial pivoting to find x from (A x = b)
/// @tparam para_dim
/// @param A May be modified inplace
/// @param b May be modified inplace
/// @param skipmask May be modified inplace
/// @param x
template<std::size_t para_dim>
inline void
GaussWithPivot(std::array<std::array<double, para_dim>, para_dim>& A,
               std::array<double, para_dim>& b,
               std::array<int, para_dim>& skipmask,
               std::array<double, para_dim>& x) {
  std::size_t maxrow;
  double maxval, maybemax;
  std::array<int, para_dim> indexmap;
  std::iota(indexmap.begin(), indexmap.end(), 0);

  // partial pivoting and forward reduction
  for (std::size_t i{0}; i < para_dim; i++) {
    // sneak in x initialization here.
    // done first, so that skipping it will have no effect on update later
    x[i] = 0.;

    // ignore clipped entries
    if (skipmask[i] != 0)
      continue;

    /* swap */
    // go through the rows and compare magnitude and mark it if needed
    maxrow = i;
    maxval = std::abs(A[i][i]);
    for (std::size_t j{i + 1}; j < para_dim; j++) {
      maybemax = std::abs(A[j][i]);
      if (maybemax > maxval) {
        maxval = maybemax;
        maxrow = j;
      }
    }
    // swap if needed. max entry row goes to i-th row.
    if (maxrow != i) {
      // swap A's rows, entries of b, indexmap, and skipmask
      std::swap(A[maxrow], A[i]);
      std::swap(b[maxrow], b[i]);
      std::swap(indexmap[maxrow], indexmap[i]);
      std::swap(skipmask[maxrow], skipmask[i]);
    }
    /* END swap */

    /* forward reduction */
    const double Aii = A[i][i];
    for (std::size_t j{i + 1}; j < para_dim; j++) {
      const double redfac = A[j][i] / Aii;
      // matrix reduction
      for (std::size_t k{i + 1}; k < para_dim; k++) {
        A[j][k] -= A[i][k] * redfac;
      }
      A[j][i] = 0.;

      // row-vec reduction0
      b[j] -= b[i] * redfac;
    }
    /* END forward reduction */
  }

  // back substitution
  double sum;
  // Looping backwards, stopping once overflow is reached, which corresponds to
  // (-1), i.e., such that the last value is 0
  for (std::size_t i{para_dim - 1}; i != static_cast<std::size_t>(-1); --i) {
    // ignore clipped entries
    if (skipmask[i] != 0)
      continue;

    sum = 0.;
    for (std::size_t j{i + 1}; j < para_dim; j++) {
      sum += A[i][j] * x[j];
    }
    x[i] = (b[i] - sum) / A[i][i];
  }

  // reorder
  CopyReorder(x, indexmap);
  CopyReorder(skipmask, indexmap);
}

} /* namespace splinepy::utils */
