#ifndef GRASS_YOSHIDA_H
#define GRASS_YOSHIDA_H

// Yoshida fourth-order integrator.

#include <complex>

namespace dyn {

template <typename F> struct Yoshida {
  std::complex<F> y0, y1;
  constexpr Yoshida() = default;
  constexpr Yoshida(std::complex<F> y0, std::complex<F> y1) : y0{y0}, y1{y1} {}

  template <typename A> void step(F h, A y2) {
    // Algorithm and constants are due to Wikipedia.

    auto constexpr CBRT2 = F(1.259921049894873164767L);
    auto constexpr W0 = -CBRT2 / (F(2) - CBRT2);
    auto constexpr W1 = F(1) / (F(2) - CBRT2);
    auto constexpr C1 = W1 / F(2), C4 = C1;
    auto constexpr C2 = (W0 + W1) / F(2), C3 = C2;
    auto constexpr D1 = W1, D3 = D1;
    auto constexpr D2 = W0;

    y0 += C1 * h * y1;
    y1 += D1 * h * y2(y0);
    y0 += C2 * h * y1;
    y1 += D2 * h * y2(y0);
    y0 += C3 * h * y1;
    y1 += D3 * h * y2(y0);
    y0 += C4 * h * y1;

    return *this;
  }
};

}; // namespace dyn

#endif // GRASS_YOSHIDA_H
