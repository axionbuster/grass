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

namespace intersect {

/// @brief Decide whether at least one intersection (point) exists between the
/// area of the disk centered at the origin with the given radius and the area
/// of the given rectangle (degenerate cases are unspecified).
/// @param ll The less-less corner of the rectangle.
/// @param gg The greater-greater corner of the rectangle.
template <typename F = float>
constexpr bool origin_disk_rectangle(F radius, std::complex<F> ll,
                                     std::complex<F> gg) noexcept {
  // https://www.jeffreythompson.org/collision-detection/circle-rect.php

  std::complex<F> test;

  if (ll.real() > F{})
    test.real(ll.real());
  else if (gg.real() < F{})
    test.real(gg.real());

  if (ll.imag() > F{})
    test.imag(ll.imag());
  else if (gg.imag() < F{})
    test.imag(gg.imag());

  return std::norm(test) < radius * radius;
}

/// @brief Decide whether at least one intersection (point) exists between a
/// circular disk and the area of a rectangle (degenerate cases are
/// unspecified).
template <typename F = float>
constexpr bool disk_rectangle(Circle<F> circ, std::complex<F> ll,
                              std::complex<F> gg) noexcept {
  // Translate the coordinate system so that the circle (circ) is at the origin.
  return origin_disk_rectangle(circ.radius, ll - circ, gg - circ);
}

} // namespace intersect

} // namespace dyn

#endif // GRASS_CIRCLE_H
