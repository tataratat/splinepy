/* Copyright (c) 2018–2021 SplineLib

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit
persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. */

template<int size>
Index<size>::Index(Length_ length, Value_ value, bool invalid) {
  using std::move;

#ifndef NDEBUG
  try {
    ThrowIfValueIsInvalid(length, value);
  } catch (OutOfRange const &exception) { Throw(exception, "splinelib::sources::utilities::Index::Index"); }
#endif

  length_ = move(length);
  value_ = move(value);
  invalid_ = move(invalid);
}

template<int size>
Index<size> & Index<size>::operator+=(Value_ const &value) {
  std_container_operations::AddAndAssignToFirst(value_, value);
#ifndef NDEBUG
  Message const kName{"splinelib::sources::utilities::Index::operator+="};

  try {
    ThrowIfValueIsInvalid(length_, value_);
  } catch (DomainError const &exception) { Throw(exception, kName); }
    catch (OutOfRange const &exception) { Throw(exception, kName); }
#endif
  return *this;
}

template<int size>
Index<size> & Index<size>::operator-=(Value_ const &value) {
  std_container_operations::SubtractAndAssignToFirst(value_, value);
#ifndef NDEBUG
  Message const kName{"splinelib::sources::utilities::Index::operator-="};

  try {
    ThrowIfValueIsInvalid(length_, value_);
  } catch (DomainError const &exception) { Throw(exception, kName); }
    catch (OutOfRange const &exception) { Throw(exception, kName); }
#endif
  return *this;
}

template<int size>
Index<size> & Index<size>::operator++() {
  for (Dimension::Type_ dimension{}; dimension < size; ++dimension) {
    Length::Type_ const &length = length_[dimension].Get();
    Index_::Type_ &value = value_[dimension].Get();
    if (length == 0) {
      continue;
    } else if (value == (length - 1)) {
      value = 0;
    } else {
      ++value;
      break;
    }
  }
  if (value_ == Value_{}) invalid_ = true;  // Invalidate the index if it was incremented from the maximum index.
  return *this;
}

template<int size>
Index<size> & Index<size>::Increment(Dimension const &dimension) {
  Dimension::Type_ const &dimension_value = dimension.Get();
#ifndef NDEBUG
  Message const kName{"splinelib::sources::utilities::Index::Increment"};

  try {
    Dimension::ThrowIfNamedIntegerIsOutOfBounds(dimension, size - 1);
  } catch (DomainError const &exception) { Throw(exception, kName, dimension_value); }
    catch (OutOfRange const &exception) { Throw(exception, kName, dimension_value); }
#endif
  Length::Type_ const &length = length_[dimension_value].Get();
  Index_::Type_ &value = value_[dimension_value].Get();
  if (length != 0) {
    if (value == (length - 1)) {
      value = 0;
    } else {
      ++value;
    }
  }
  return *this;
}

template<int size>
Index<size> & Index<size>::operator--() {
  if (value_ == Value_{}) {  // Invalidate the index if it is decremented from the zero index.
    invalid_ = true;
  } else {
    for (Dimension::Type_ dimension{}; dimension < size; ++dimension) {
      Length::Type_ const &length = length_[dimension].Get();
      Index_::Type_ &value = value_[dimension].Get();
      if (length == 0) {
        continue;
      } else if (value == 0) {
        value = (length - 1);
      } else {
        --value;
        break;
      }
    }
  }
  return *this;
}

template<int size>
Index<size> & Index<size>::Decrement(Dimension const &dimension) {
  Dimension::Type_ const &dimension_value = dimension.Get();
#ifndef NDEBUG
  Message const kName{"splinelib::sources::utilities::Index::Decrement"};

  try {
    Dimension::ThrowIfNamedIntegerIsOutOfBounds(dimension, size - 1);
  } catch (DomainError const &exception) { Throw(exception, kName, dimension_value); }
    catch (OutOfRange const &exception) { Throw(exception, kName, dimension_value); }
#endif
  Length::Type_ const &length = length_[dimension_value].Get();
  Index_::Type_ &value = value_[dimension_value].Get();
  if (length != 0) {
    if (value == 0) {
      value = (length - 1);
    } else {
      --value;
    }
  }
  return *this;
}

template<int size>
Index<size> operator+(Index<size> const &lhs, typename Index<size>::Value_ const &rhs) {
  Index sum{lhs};
#ifndef NDEBUG
  Message const kName{"splinelib::sources::utilities::operator+::Index"};

  try {
#endif
    sum += rhs;
#ifndef NDEBUG
  } catch (DomainError const &exception) { Throw(exception, kName); }
    catch (OutOfRange const &exception) { Throw(exception, kName); }
#endif
  return sum;
}

template<int size>
Index<size> operator-(Index<size> const &lhs, typename Index<size>::Value_ const &rhs) {
  Index difference{lhs};
#ifndef NDEBUG
  Message const kName{"splinelib::sources::utilities::operator-::Index"};

  try {
#endif
    difference -= rhs;
#ifndef NDEBUG
  } catch (DomainError const &exception) { Throw(exception, kName); }
    catch (OutOfRange const &exception) { Throw(exception, kName); }
#endif
  return difference;
}

template<int size>
bool operator==(Index<size> const &lhs, Index<size> const &rhs) {
  return ((lhs.invalid_ == rhs.invalid_) && (lhs.length_ == rhs.length_) && (lhs.value_ == rhs.value_));
}

template<int size>
bool operator!=(Index<size> const &lhs, Index<size> const &rhs) {
  return (!(lhs == rhs));
}

template<int size>
splinelib::Index const & Index<size>::operator[](Dimension const &dimension) const {
  Dimension::Type_ const &dimension_value = dimension.Get();
#ifndef NDEBUG
  Message const kName{"splinelib::sources::utilities::Index::operator[]"};

  try {
    Dimension::ThrowIfNamedIntegerIsOutOfBounds(dimension, size - 1);
  } catch (DomainError const &exception) { Throw(exception, kName, dimension_value); }
    catch (OutOfRange const &exception) { Throw(exception, kName, dimension_value); }
#endif
  return value_[dimension_value];
}

template<int size>
Index<size> Index<size>::First(Length_ length) {
  return Index{std::move(length)};
}

template<int size>
Index<size> Index<size>::Last(Length_ length) {
  using std::move;

  Value_ value;
  std::transform(length.begin(), length.end(), value.begin(),
                 [] (Length const &length_for_dimension) { Index_::Type_ const &length = length_for_dimension.Get();
                     if (length == 0) { return Index_{}; } else { return Index_{length - 1}; } });
  return Index{move(length), move(value), };
}

template<int size>
Index<size> Index<size>::Behind(Length_ length) {
  return ++Last(std::move(length));
}

template<int size>
Index<size> Index<size>::Before(Length_ length) {
  return --First(std::move(length));
}

template<int size>
int Index<size>::GetTotalNumberOfIndices() const {
  return DetermineStride(length_, Dimension{size});
}

template<int size>
typename Index<size>::Value_ Index<size>::GetIndex() const {
  return value_;
}

template<int size>
splinelib::Index Index<size>::GetIndex1d() const {
  Index_ index_1d{};
  Dimension::ForEach(0, size, [&] (Dimension const &dimension) {
      index_1d += (value_[dimension.Get()] * DetermineStride(length_, dimension)); });
  return index_1d;
}

template<int size>
int Index<size>::DetermineStride(Length_ const &length, Dimension const &dimension) const {
  typename Length_::const_iterator const &begin = length.begin();
  return std::transform_reduce(begin, begin + dimension.Get(), 1, std::multiplies{}, [] (Length const &length) {
                                   if (length == Length{}) { return 1; } else { return length.Get(); } });
}

#ifndef NDEBUG
template<int size>
void Index<size>::ThrowIfValueIsInvalid(Length_ const &length, Value_ const &value) {
  Dimension::ForEach(0, size, [&] (Dimension const &dimension) {
      Dimension::Type_ const &current_dimension = dimension.Get();
      Message const dimension_string{"for dimension " + std::to_string(current_dimension) + ": "};
      try {
        Index_::ThrowIfNamedIntegerIsOutOfBounds(value[current_dimension],
                                                 std::max(0, length[current_dimension].Get() - 1));
      } catch (DomainError const &exception) { throw DomainError(dimension_string + exception.what()); }
        catch (OutOfRange const &exception) { throw OutOfRange(dimension_string + exception.what()); }});
}
#endif
