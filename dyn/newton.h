#ifndef GRASS_NEWTON_H
#define GRASS_NEWTON_H

#include <array>
#include <cmath>
#include <complex>

#include "circle.h"
#include "halton.h"
#include "kahan.h"

namespace dyn {

template <typename F = float, unsigned short N_MONTE = 30> class Gravity {
  std::array<std::complex<F>, 2 * N_MONTE> disk{};
  F G{1};

public:
  Gravity() { init_disk(); }
  explicit Gravity(F G) : G{G} { init_disk(); }

  std::complex<F> field(Circle<F> c0, Circle<F> c1, F m1) noexcept {
    c1 -= c0, c0 -= c0;
    auto s = c1.radius + c0.radius, d = std::abs(c1.radius - c0.radius),
         r = std::abs(c1), t = F(1) / r;
    if (s <= r)
      // Disjoint
      return t * t * t * m1 * G * c1;
    else if (d <= r)
      // Intersecting
      return when_intersecting(c0.radius, c1, m1);
    else
      // Contains (no force)
      return 0;
  }

private:
  constexpr void init_disk() {
    Halton<F, 2> h2;
    Halton<F, 3> h3;
    for (unsigned short i = 0; i < N_MONTE; i++)
      disk[i] = dyn::disk(std::complex<F>{h2.xy(), h3.xy()});
    for (unsigned short i = 0; i < N_MONTE; i++)
      disk[N_MONTE + i] = std::conj(disk[i]);
  }

  std::complex<F> when_intersecting(F r0, Circle<F> c1, F m1) noexcept {
    Kahan<F> B;
    for (auto &&p : disk) {
      auto q = c1 - p;
      auto r = std::abs(q);
      if (r <= c1.radius)
        continue;
      auto s = F(1) / r;
      auto dB = s * s * s * q;
      B += dB;
    }
    auto constexpr premultiply = F(1) / F(N_MONTE);
    return premultiply * m1 * B();
  }
};

} // namespace dyn

#endif // GRASS_NEWTON_H
