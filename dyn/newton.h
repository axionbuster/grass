#ifndef GRASS_NEWTON_H
#define GRASS_NEWTON_H

/// @file newton.h
/// @brief Computation of Newtonian gravity.

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>

#include "circle.h"
#include "halton.h"

namespace dyn {

/// @brief Compute the Newtonian gravitational interaction between pairs of
/// circles, taking into account when they are too close to one another.
/// @tparam F A floating-point type.
/// @tparam N_MONTE Number of Monte Carlo trials.
template <typename F = float, unsigned short N_MONTE = 30> class Gravity {
  /// @brief Quasi-random points on the unit disk centered about the origin used
  /// for Monte Carlo integration in the case of overlapping circles.
  std::array<std::complex<F>, N_MONTE> disk{};

  // Each Halton sequence (a kind of low-discrepancy sequence) creates an
  // evenly spaced set of points on the unit interval (0, 1); unlike the
  // uniform distribution, however, the points look "uniformly distributed"
  // (number of points being mostly proportional to length of any subset)
  // even for a finite sample of points.
  Halton<F, 2> h2;
  Halton<F, 3> h3;

public:
  /// @brief Create an instance with a quasi-random internal state.
  constexpr Gravity() { refresh_disk(); }

  /// @brief Compute the gravitational attraction that a test particle
  /// represented by the circle c0 due to a mass of circle c1 and mass m1.
  /// @param distance Optional distance (or zero if needs to be computed).
  /// @return A vector quantity of dimensions [M/L/L]. Divide it by the
  /// universal gravitational constant ("G") to get L/T/T.
  std::complex<F> field(Circle<F> c0, Circle<F> c1, F m1,
                        F distance = F{}) const noexcept {
    // Use complex arithmetic to translate the coordinate system so that c0
    // appears to be at the origin, but let the respective radii be unaffected.
    c1 -= c0;

    // Are the circles...
    auto s = c1.radius + c0.radius, d = std::abs(c1.radius - c0.radius),
         r = distance ? distance : std::abs(c1), t = F(1) / r;
    if (s <= r)
      // Disjoint? Use the usual law, treating these circles as point
      // particles whose masses are concentrated at the given centers.
      return t * t * t * m1 * c1;
    else if (d <= r)
      // Intersecting?
      return when_intersecting(c0.radius, c1, m1);
    else
      // Or, does one circle fully contain the other? Then, apply the shell
      // theorem to get a net force of zero.
      return F(0);
  }

  /// @brief Populate internal random disk (used for calculating forces in the
  /// case of intersecting circles) with new evenly distributed points on the
  /// unit disk centered about the origin. Call often to avoid bias.
  constexpr void refresh_disk() {
    // Fill `disk` with random points on unit disk centered about origin
    // by rejection sampling.
    for (auto &&p : disk) {
      do {
        // Scale and move (0,1) x (0,1) square to (-1,1) x (-1,1) square.
        p = F(2) * std::complex<F>{h2.x01(), h3.x01()} -
            std::complex<F>{F(1), F(1)};
        // Use of `norm` makes this function constexpr (where `abs` or
        // trigonometric functions are not allowed as of C++20).
      } while (std::norm(p) >= F(1));
    }

    // Attempt to improve branch prediction somewhat by sorting the points about
    // some axis (here, the real axis).
    std::sort(disk.begin(), disk.end(),
              [](auto p, auto q) { return p.real() < q.real(); });
  }

private:
  /// @brief Compute the gravitational attraction that a test particle (center
  /// of radius r0) at the origin feels due to the presence of source particle
  /// (represented by the circle c1) with mass m1.
  /// @param r0 Radius of the test particle (at the origin).
  /// @param c1 Center and radius of source particle.
  /// @param m1 Mass of particle c1.
  /// @return Vector quantity of dimensions M/L/L.
  std::complex<F> when_intersecting(F r0, Circle<F> c1, F m1) const noexcept {
    // Assuming uniform mass distribution (by area), divide the test particle's
    // circle into many small pieces. To each particle, apply Newton's shell
    // theorem: (a) Inside a radially symmetrical body, no force is felt due to
    // that body; (b) outside it, the force is as though its mass [m1] was
    // concentrated at the center of it [c1].

    // Unweighted inverse-square summation (B).
    std::complex<F> B;
    for (auto &&p : disk) {
      auto q = c1 - r0 * p;
      auto r = std::abs(q);
      if (r > c1.radius) {
        auto s = F(1) / r;
        B += s * s * s * q;
      }
    }

    // Converting constexpr integer N_MONTE to floating point constexpr is
    // helpful because it reduces burden due to type conversion on CPU.
    return F(1) / F(N_MONTE) * m1 * B;
  }
};

} // namespace dyn

#endif // GRASS_NEWTON_H
