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

/// @brief Decide whether at least one intersection (point) exists between a
/// circular disk and a rectangle (degenerate cases are unspecified).
template <typename F = float>
constexpr bool disk_rect_intersect(Circle<F> circ, std::complex<F> ll,
                                   std::complex<F> gg) {
  ll -= circ, gg -= circ;
  auto rad = circ.radius;

  // Suppose that the circle and the rectangle intersect, then:

  // 1. The circle's bounding square and the rectangle must intersect.
  if (ll.real() > rad || ll.imag() > rad || gg.real() < -rad ||
      gg.imag() < -rad)
    // (They don't. -> No intersection.)
    return false;
  // 2. Otherwise, if either side of the rectangle cut through a coordinate
  // axis then they intersect.
  if (std::signbit(ll.real()) != std::signbit(gg.real()) ||
      std::signbit(ll.imag()) != std::signbit(gg.imag()))
    // (It does. -> Intersection.)
    return true;
  // 3. Lastly, otherwise, if any of the four corners are in the disk bounded
  // by the circle then the circle and the rectangle intersect.
  F cc[4] = {std::abs(ll), std::abs(gg),
             std::abs(std::complex<F>{ll.real(), gg.imag()}),
             std::abs(std::complex<F>{ll.imag(), gg.real()})};
  for (auto &&c : cc)
    if (c < rad)
      return true;
  // No other intersecting cases exist.
  return false;
}

} // namespace dyn

#endif // GRASS_CIRCLE_H
