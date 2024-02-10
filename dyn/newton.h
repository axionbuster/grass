#ifndef GRASS_NEWTON_H
#define GRASS_NEWTON_H

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>

#include "circle.h"
#include "halton.h"

namespace dyn {

/// @brief Compute the Newtonian gravitational interaction between pairs of
/// circles, taking into account when they are too close to one another.
template <typename F = float, unsigned short N_MONTE = 30> class Gravity {
  /// @brief Quasi-random points on the unit disk centered about the origin used
  /// for Monte Carlo integration in the case of overlapping circles.
  std::array<std::complex<F>, N_MONTE> disk{};

public:
  /// @brief Universal gravitational constant.
  F G{1};

  /// @brief Create an instance with a quasi-random internal state.
  constexpr Gravity() { init_disk(); }

  /// @brief Compute the acceleration that a test particle represented by the
  /// circle c0 due to a mass of circle c1 and mass m1.
  std::complex<F> field(Circle<F> c0, Circle<F> c1, F m1) const noexcept {
    // Use complex arithmetic to translate the coordinate system so that c0
    // appears to be at the origin, but let the respective radii be unaffected.
    c1 -= c0, c0 -= c0;

    // Are the circles...
    auto s = c1.radius + c0.radius, d = std::abs(c1.radius - c0.radius),
         r = std::abs(c1), t = F(1) / r;
    if (s <= r)
      // Disjoint?
      return t * t * t * m1 * G * c1;
    else if (d <= r)
      // Intersecting?
      return when_intersecting(c0.radius, c1, m1);
    else
      // Or, one circle fully contains other?
      return 0;
  }

private:
  /// @brief Populate `disk` with evenly distributed points on the unit disk
  /// centered about the origin.
  constexpr void init_disk() {
    // Halton sequences make an evenly distributed grid of points.
    Halton<F, 2> h2;
    Halton<F, 3> h3;

    // Skip the first few terms.
    for (int i = 0; i < 1234; i++)
      h2.x01(), h3.x01();

    // Fill `disk` with random points on unit disk centered about origin
    // by rejection sampling.
    for (auto &&p : disk) {
      do {
        p = F(2) * std::complex<F>{h2.x01(), h3.x01()} -
            std::complex<F>{F(1), F(1)};
      } while (std::norm(p) >= F(1));
    }

    // Attempt to improve branch prediction somewhat by sorting the points about
    // some axis (here, the real axis).
    std::sort(disk.begin(), disk.end(),
              [](auto p, auto q) { return p.real() < q.real(); });
  }

  std::complex<F> when_intersecting(F r0, Circle<F> c1, F m1) const noexcept {
    // Chop "c0" (here, c0 is centered about the origin, and only the radius
    // r0 is relevant) into small pieces, assuming that the mass is proportional
    // to the area of each piece. Then, add up the small forces for the pieces
    // that are outside c1, but ignore the pieces inside c1, by application of
    // Newton's shell theorem, which says that, within a radially symmetric
    // shell of a body, all gravitational forces due to that body cancel out.

    // Unweighted inverse-square summation (B).
    std::complex<F> B;
    for (auto &&p : disk) {
      auto q = c1 - r0 * p;
      auto r = std::abs(q);
      if (r > c1.radius) {
        // Give up a bit of accuracy for speed.
        auto s = F(1) / r;
        auto dB = s * s * s * (c1 - r0 * p);
        B += dB;
      }
    }

    // Converting constexpr integer N_MONTE to floating point constexpr is
    // helpful because it reduces burden due to type conversion on CPU.
    return F(0.5) / F(N_MONTE) * m1 * G * B;
  }
};

} // namespace dyn

#endif // GRASS_NEWTON_H
