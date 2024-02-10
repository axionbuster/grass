#ifndef GRASS_CIRCLE_H
#define GRASS_CIRCLE_H

#include <complex>

namespace dyn {

template <typename F = float> struct Circle : public std::complex<F> {
  F radius{1};
  constexpr Circle() = default;
  constexpr Circle(std::complex<F> center, F radius)
      : std::complex<F>(center), radius(radius) {}
};

} // namespace dyn

#endif // GRASS_CIRCLE_H
