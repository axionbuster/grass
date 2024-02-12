#ifndef GRASS_CIRCLE_H
#define GRASS_CIRCLE_H

/// @file circle.h
/// @brief Represent a circle by the center and radius.

#include <complex>
#include <numbers>

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

/// @brief Decide whether at least one intersection (point) exists between the
/// area of the disk centered at the origin with the given radius and the area
/// of the given rectangle.
/// @param ll The less-less corner of the rectangle.
/// @param gg The greater-greater corner of the rectangle.
template <typename F = float>
bool origin_disk_rect_isct(F radius, std::complex<F> ll,
                           std::complex<F> gg) noexcept {
  // sin(45 deg) [also cos(45 deg)].
  constexpr F sc45 = std::numbers::sqrt2_v<F> / F(2);

  // The rectangle's four corners are in four different quadrants.
  if (std::signbit(ll.real()) != std::signbit(gg.real()) &&
      std::signbit(ll.imag()) != std::signbit(gg.imag()))
    return true;

  // Decide whether the given rectangle and a square centered at (0,0) having a
  // "radius" (half side length) of `radius` (thus a full side length of twice
  // the `radius`).
  auto recsqisct = [&]() {
    return -radius < gg.real() && ll.real() < radius && -radius < gg.imag() &&
           gg.imag() < radius;
  };

  // The rectangle does not touch the bounding square of the circle.
  if (!recsqisct())
    return false;

  // The rectangle touches the interior square of the circle.
  radius *= sc45;
  if (recsqisct())
    return true;

  // They don't touch.
  return false;
}

/// @brief Decide whether at least one intersection (point) exists between a
/// circular disk and a rectangle (degenerate cases are unspecified).
template <typename F = float>
bool disk_rect_intersect(Circle<F> circ, std::complex<F> ll,
                         std::complex<F> gg) noexcept {
  ll -= circ, gg -= circ;
  return origin_disk_rect_isct(circ.radius, ll, gg);
}

} // namespace dyn

#endif // GRASS_CIRCLE_H
