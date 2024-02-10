#ifndef GRASS_YOSHIDA_H
#define GRASS_YOSHIDA_H

/// @file yoshida.h
/// @brief Yoshida's fourth-order symplectic (area-preserving) integrator.

#include <complex>

namespace dyn {

/// @brief Yoshida's fourth-order symplectic (area-preserving) integrator for
/// complex numbers. "Area-preserving" integrators will preserve the energy of a
/// system of differential equations.
/// @tparam F A floating-point type.
template <typename F = float> struct Yoshida {
  /// @brief The zeroth and first derivatives, respectively.
  std::complex<F> y0, y1;

  /// @brief Create a zero-initialized instance.
  constexpr Yoshida() = default;

  /// @brief Instantiate with given zeroth and first derivative values.
  constexpr Yoshida(std::complex<F> y0, std::complex<F> y1) : y0{y0}, y1{y1} {}

  /// @brief Compute the second derivative three times with slightly different
  /// zeroth-derivative values (y0) and update both the zeroth- and
  /// first-derivative values of the internal state.
  /// @param y2 An effectively stateless function that takes in a complex
  /// zeroth-derivative value and computes the complex second-derivative.
  /// @param h Step size (finite, positive number).
  template <typename Y2> void step(F h, Y2 y2) {
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
  }
};

}; // namespace dyn

#endif // GRASS_YOSHIDA_H
