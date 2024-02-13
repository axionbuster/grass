#ifndef GRASS_VERLET_H
#define GRASS_VERLET_H

#include <complex>

namespace dyn {

/// @brief Velocity Verlet integration.
template <typename F = float> struct Verlet {
  /// @brief Last values of the zeroth and first derivatives, respectively.
  std::complex<F> y0, y1;

  /// @brief Compute the second derivative twice and then update the zeroth
  /// and first derivatives, respectively.
  /// @param h Step size.
  /// @param y2 A function that takes in a complex zeroth derivative value and
  /// then returns a complex second derivative.
  template <typename A> void step(F h, A y2) {
    auto a = y2(y0);
    y0 += h * y1 + h * h * F(0.5) * a;
    auto b = y2(y0);
    y1 += h * F(0.5) * (a + b);
  }
};

} // namespace dyn

#endif // GRASS_VERLET_H
