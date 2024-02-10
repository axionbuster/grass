#ifndef GRASS_KAHAN_H
#define GRASS_KAHAN_H

/// @file kahan.h
/// @brief Kahan's compensated summation.

namespace dyn {

/// @brief Kahan's compensated summation. Subnormal numbers must be enabled.
/// @tparam T A type that supports initialization, addition, subtraction, and
/// copy.
template <typename T = float> class Kahan {
  /// @brief The accumulator.
  T a{};

  /// @brief The error.
  T e{};

public:
  /// @brief Construct a zero-initialized instance.
  Kahan() = default;

  /// @brief Construct an instance with a value for the accumulator.
  /// @param a Initial value of the accumulator.
  explicit Kahan(T a) : a{a} {};

  /// @brief Return the accumulator.
  /// @return The accumulator.
  constexpr T operator()() const { return a; }

  /// @brief Add a value to the accumulator and update the error term.
  constexpr Kahan &add(T v) {
    // Precondition:
    // - For any two T values t and s and a certain T-type constant "0," it must
    // be that (t-s == 0) if and only if (t == s).
    // - If T is a floating-point type, this condition implies the existence of
    // subnormal numbers.
    auto y = v - e, t = a + y;
    return e = t - a - y, a = t, *this;
  }

  /// @brief Add a value to the accumulator and update the error term.
  constexpr Kahan &operator+=(T v) { return add(v); }
};

} // namespace dyn

#endif // GRASS_KAHAN_H
