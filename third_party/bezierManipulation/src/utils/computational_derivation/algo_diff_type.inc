#ifndef UTILS_COMPUTATIONAL_DERIVATION_ALGO_DIFF_TYPE_INC
#define UTILS_COMPUTATIONAL_DERIVATION_ALGO_DIFF_TYPE_INC

// Addition
template <typename Scalar>
constexpr AlgoDiffType<Scalar> AlgoDiffType<Scalar>::operator+(
    const AlgoDiffType<Scalar>& b) const {
  ADT_ result_value{(*this)};
  result_value += b;
  return result_value;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar>& AlgoDiffType<Scalar>::operator+=(
    const AlgoDiffType<Scalar>& b) {
  assert(b.GetNumberOfDerivatives() == GetNumberOfDerivatives());
  v_ += b.v_;
  for (IndexingType_ i{}; i < GetNumberOfDerivatives(); i++) {
    d_[i] += b.d_[i];
  }
  return (*this);
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> AlgoDiffType<Scalar>::operator+(
    const Scalar& b) const {
  ADT_ result_value{(*this)};
  result_value += b;
  return result_value;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar>& AlgoDiffType<Scalar>::operator+=(
    const Scalar& b) {
  v_ += b;
  return (*this);
}

// Substraction and negation
template <typename Scalar>
constexpr AlgoDiffType<Scalar> AlgoDiffType<Scalar>::operator-() const {
  ADT_ result_value{(*this)};
  result_value.v_ = -result_value.v_;
  for (IndexingType_ i{}; i < GetNumberOfDerivatives(); i++) {
    result_value.d_[i] = -result_value.d_[i];
  }
  return result_value;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> AlgoDiffType<Scalar>::operator-(
    const AlgoDiffType<Scalar>& b) const {
  assert(b.GetNumberOfDerivatives() == GetNumberOfDerivatives());
  ADT_ result_value{(*this)};
  result_value -= b;
  return result_value;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar>& AlgoDiffType<Scalar>::operator-=(
    const AlgoDiffType<Scalar>& b) {
  v_ -= b.v_;
  for (IndexingType_ i{}; i < GetNumberOfDerivatives(); i++) {
    d_[i] -= b.d_[i];
  }
  return (*this);
}
template <typename Scalar>
constexpr AlgoDiffType<Scalar> AlgoDiffType<Scalar>::operator-(
    const Scalar& b) const {
  ADT_ result_value{(*this)};
  result_value -= b;
  return result_value;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar>& AlgoDiffType<Scalar>::operator-=(
    const Scalar& b) {
  v_ -= b;
  return (*this);
}

// Multiplication
template <typename Scalar>
constexpr AlgoDiffType<Scalar> AlgoDiffType<Scalar>::operator*(
    const AlgoDiffType<Scalar>& b) const {
  ADT_ result_value{(*this)};
  result_value *= b;
  return result_value;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar>& AlgoDiffType<Scalar>::operator*=(
    const AlgoDiffType<Scalar>& b) {
  assert(b.GetNumberOfDerivatives() == GetNumberOfDerivatives());
  for (IndexingType_ i{}; i < GetNumberOfDerivatives(); i++) {
    d_[i] = b.d_[i] * v_ + b.v_ * d_[i];  // f'(x) = u'*v +  u*v'
  }
  // Must be at the end otherwise overwrite value
  v_ *= b.v_;
  return (*this);
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> AlgoDiffType<Scalar>::operator*(
    const Scalar& b) const {
  ADT_ result_value{(*this)};
  result_value *= b;
  return result_value;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar>& AlgoDiffType<Scalar>::operator*=(
    const Scalar& b) {
  v_ *= b;
  for (IndexingType_ i{}; i < GetNumberOfDerivatives(); i++) {
    d_[i] *= b;
  }
  return (*this);
}

// Division
template <typename Scalar>
constexpr AlgoDiffType<Scalar> AlgoDiffType<Scalar>::operator/(
    const AlgoDiffType<Scalar>& b) const {
  ADT_ result_value{(*this)};
  result_value /= b;
  return result_value;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar>& AlgoDiffType<Scalar>::operator/=(
    const AlgoDiffType<Scalar>& b) {
  assert(b.GetNumberOfDerivatives() == GetNumberOfDerivatives());
  const Scalar_ inverse_b_squared = 1.0 / (b.v_ * b.v_);
  for (IndexingType_ i{}; i < GetNumberOfDerivatives(); i++) {
    // f'(x) = (u'*v -  u*v') / v^2
    d_[i] = (b.v_ * d_[i] - b.d_[i] * v_) * inverse_b_squared;
  }
  v_ /= b.v_;
  return (*this);
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> AlgoDiffType<Scalar>::operator/(
    const Scalar& b) const {
  ADT_ result_value{(*this)};
  result_value /= b;
  return result_value;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar>& AlgoDiffType<Scalar>::operator/=(
    const Scalar& b) {
  const Scalar inverse_b{1 / b};
  v_ *= inverse_b;
  for (IndexingType_ i{}; i < GetNumberOfDerivatives(); i++) {
    d_[i] *= inverse_b;
  }
  return (*this);
}

//////////////////////////////
// Friend Functions
//////////////////////////////

template <typename Scalar>
std::ostream& operator<<(std::ostream& os, const AlgoDiffType<Scalar>& a) {
  os << a.v_;
  return os;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> operator+(const Scalar& a,
                                         const AlgoDiffType<Scalar>& b) {
  return b + a;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> operator-(const Scalar& a,
                                         const AlgoDiffType<Scalar>& b) {
  return (-b) + a;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> operator*(const Scalar& a,
                                         const AlgoDiffType<Scalar>& b) {
  return b * a;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> operator/(const Scalar& a,
                                         const AlgoDiffType<Scalar>& b) {
  const AlgoDiffType<Scalar> one_type{a, b.GetNumberOfDerivatives()};
  // A few too many operations (multiplication with 0), but typesafe
  return one_type / b;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> exp(const AlgoDiffType<Scalar>& exponent) {
  Scalar tmp;
  if constexpr (std::is_arithmetic_v<Scalar>) {
    // Use STD namespace function
    tmp = std::exp(exponent.v_);
  } else {
    tmp = exp(exponent.v_);
  }
  AlgoDiffType<Scalar> result_value{exponent};
  result_value.v_ = Scalar{1.};
  result_value *= tmp;
  return result_value;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> abs(const AlgoDiffType<Scalar>& base) {
  return base.v_ > 0.0 ? base : (-base);
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> pow(const AlgoDiffType<Scalar>& base,
                                   const Scalar& power) {
  Scalar tmp;
  if constexpr (std::is_arithmetic_v<Scalar>) {
    // Use STD namespace function
    tmp = std::pow(base.v_, power - 1.0);
  } else {
    tmp = pow(base.v_, power - 1.0);
  }
  AlgoDiffType<Scalar> result_value{base};
  result_value.v_ *= tmp;
  for (std::size_t i{}; i < base.GetNumberOfDerivatives(); i++) {
    result_value.d_[i] *= power * tmp;
  }
  return result_value;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> pow(const AlgoDiffType<Scalar>& base,
                                   const AlgoDiffType<Scalar>& power) {
  return exp(log(base) * power);
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> sqrt(const AlgoDiffType<Scalar>& radicand) {
  Scalar root;
  if constexpr (std::is_arithmetic_v<Scalar>) {
    // Use STD namespace function
    root = std::sqrt(radicand.v_);
  } else {
    root = sqrt(radicand.v_);
  }
  const Scalar& half_inverse_root = 0.5 / root;
  AlgoDiffType<Scalar> result_value{radicand};
  result_value.v_ = root;
  for (std::size_t i{}; i < radicand.GetNumberOfDerivatives(); i++) {
    result_value.d_[i] *= half_inverse_root;
  }
  return result_value;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> log(const AlgoDiffType<Scalar>& xi) {
  using namespace std;               // required for intrinsic types
  const Scalar& temp = 1.0 / xi.v_;  // Multiplication is cheaper than division
  AlgoDiffType<Scalar> result_value{xi};
  if constexpr (std::is_arithmetic_v<Scalar>) {
    // Use STD namespace function
    result_value.v_ = std::log(xi.v_);
  } else {
    result_value.v_ = log(xi.v_);
  }
  for (std::size_t i{}; i < xi.GetNumberOfDerivatives(); i++) {
    result_value.d_[i] *= temp;
  }
  return result_value;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> log10(const AlgoDiffType<Scalar>& a) {
  using namespace std;  // requ ired for intrinsic types
  const Scalar& temp = 1.0 / (a.v_ * std::log(10.));
  AlgoDiffType<Scalar> result_value{a};
  if constexpr (std::is_arithmetic_v<Scalar>) {
    // Use STD namespace function
    result_value.v_ = std::log10(a.v_);
  } else {
    result_value.v_ = log10(a.v_);
  }
  for (std::size_t i{}; i < a.GetNumberOfDerivatives(); i++) {
    result_value.d_[i] *= temp;
  }
  return result_value;
};

/////////
// Trigonometric functions

template <typename Scalar>
constexpr AlgoDiffType<Scalar> cos(const AlgoDiffType<Scalar>& angle) {
  AlgoDiffType<Scalar> result_value{angle};
  Scalar sin_of_angle;
  if constexpr (std::is_arithmetic_v<Scalar>) {
    // Use STD namespace function
    sin_of_angle = -std::sin(angle.v_);
    result_value.v_ = std::cos(angle.v_);
  } else {
    sin_of_angle = -sin(angle.v_);
    result_value.v_ = cos(angle.v_);
  }
  for (std::size_t i{}; i < angle.GetNumberOfDerivatives(); i++) {
    result_value.d_[i] *= sin_of_angle;
  }
  return result_value;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> sin(const AlgoDiffType<Scalar>& angle) {
  AlgoDiffType<Scalar> result_value{angle};
  Scalar cos_of_angle;
  if constexpr (std::is_arithmetic_v<Scalar>) {
    // Use STD namespace function
    cos_of_angle = std::cos(angle.v_);
    result_value.v_ = std::sin(angle.v_);
  } else {
    cos_of_angle = cos(angle.v_);
    result_value.v_ = sin(angle.v_);
  }
  for (std::size_t i{}; i < angle.GetNumberOfDerivatives(); i++) {
    result_value.d_[i] *= cos_of_angle;
  }
  return result_value;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> tan(const AlgoDiffType<Scalar>& angle) {
  AlgoDiffType<Scalar> result_value{angle};
  Scalar auxiliary_value;
  if constexpr (std::is_arithmetic_v<Scalar>) {
    // Use STD namespace function
    auxiliary_value = 1. / std::cos(angle.v_);
    auxiliary_value *= auxiliary_value;
    result_value.v_ = std::tan(angle.v_);
  } else {
    auxiliary_value = 1. / cos(angle.v_);
    auxiliary_value *= auxiliary_value;
    result_value.v_ = tan(angle.v_);
  }
  for (std::size_t i{}; i < angle.GetNumberOfDerivatives(); i++) {
    result_value.d_[i] *= auxiliary_value;
  }
  return result_value;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> acos(const AlgoDiffType<Scalar>& angle) {
  AlgoDiffType<Scalar> result_value{angle};
  Scalar auxiliary_value;
  if constexpr (std::is_arithmetic_v<Scalar>) {
    // Use STD namespace function
    auxiliary_value = -1. / std::sqrt(1 - angle.v_ * angle.v_);
    result_value.v_ = std::acos(angle.v_);
  } else {
    auxiliary_value = -1. / sqrt(1 - angle.v_ * angle.v_);
    result_value.v_ = cos(angle.v_);
  }
  for (std::size_t i{}; i < angle.GetNumberOfDerivatives(); i++) {
    result_value.d_[i] *= auxiliary_value;
  }
  return result_value;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> asin(const AlgoDiffType<Scalar>& angle) {
  AlgoDiffType<Scalar> result_value{angle};
  Scalar auxiliary_value;
  if constexpr (std::is_arithmetic_v<Scalar>) {
    // Use STD namespace function
    auxiliary_value = 1. / std::sqrt(1 - angle.v_ * angle.v_);
    result_value.v_ = std::asin(angle.v_);
  } else {
    auxiliary_value = 1. / std::sqrt(1 - angle.v_ * angle.v_);
    result_value.v_ = asin(angle.v_);
  }
  for (std::size_t i{}; i < angle.GetNumberOfDerivatives(); i++) {
    result_value.d_[i] *= auxiliary_value;
  }
  return result_value;
}

template <typename Scalar>
constexpr AlgoDiffType<Scalar> atan(const AlgoDiffType<Scalar>& angle) {
  AlgoDiffType<Scalar> result_value{angle};
  const Scalar auxiliary_value{1. / (1. + angle.v_ * angle.v_)};
  if constexpr (std::is_arithmetic_v<Scalar>) {
    // Use STD namespace function
    result_value.v_ = std::atan(angle.v_);
  } else {
    result_value.v_ = atan(angle.v_);
  }
  for (std::size_t i{}; i < angle.GetNumberOfDerivatives(); i++) {
    result_value.d_[i] *= auxiliary_value;
  }
  return result_value;
}

template <typename Scalar>
constexpr bool operator>(const Scalar& scalar,
                         const AlgoDiffType<Scalar>& adt) {
  return adt < scalar;
};

template <typename Scalar>
constexpr bool operator>=(const Scalar& scalar,
                          const AlgoDiffType<Scalar>& adt) {
  return adt <= scalar;
};

template <typename Scalar>
constexpr bool operator<(const Scalar& scalar,
                         const AlgoDiffType<Scalar>& adt) {
  return adt > scalar;
};

template <typename Scalar>
constexpr bool operator<=(const Scalar& scalar,
                          const AlgoDiffType<Scalar>& adt) {
  return adt >= scalar;
};

template <typename Scalar>
constexpr bool operator==(const Scalar& scalar,
                          const AlgoDiffType<Scalar>& adt) {
  return adt == scalar;
};

template <typename Scalar>
constexpr bool operator!=(const Scalar& scalar,
                          const AlgoDiffType<Scalar>& adt) {
  return !(adt == scalar);
};

#endif  // UTILS_COMPUTATIONAL_DERIVATION_ALGO_DIFF_TYPE_INC
