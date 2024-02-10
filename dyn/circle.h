#ifndef GRASS_CIRCLE_H
#define GRASS_CIRCLE_H

/// @file circle.h
/// @brief Represent a circle by the center and radius.

#include <complex>

namespace dyn {

/// @brief A circle (center and radius). Complex arithmetic is applied to
/// centers but not to radii.
/// @tparam F A floating-point type.
template <typename F = float> struct Circle : public std::complex<F> {
  F radius{1};

  /// @brief Construct a unit circle about the origin.
  constexpr Circle() = default;

  /// @brief Construct a circle with the given center and radius.
  constexpr Circle(std::complex<F> center, F radius)
      : std::complex<F>(center), radius(radius) {}
};

} // namespace dyn

#endif // GRASS_CIRCLE_H
