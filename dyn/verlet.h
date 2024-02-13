#ifndef GRASS_VERLET_H
#define GRASS_VERLET_H

#include <complex>

namespace dyn {

template <typename F = float> struct Verlet {
  std::complex<F> y0, y1;

  template <typename A> void step(F h, A y2) {
    auto a = y2(y0);
    y0 += h * y1 + h * h * F(0.5) * a;
    auto b = y2(y0);
    y1 += h * F(0.5) * (a + b);
  }
};

} // namespace dyn

#endif // GRASS_VERLET_H
