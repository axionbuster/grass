#ifndef GRASS_HALTON_H
#define GRASS_HALTON_H

#include <complex>
#include <numbers>

namespace dyn {

template <typename F = float, unsigned short b = 2, unsigned int LIM = 0x1000>
class Halton {
  unsigned short i{};

public:
  constexpr static F xy(unsigned short i) {
    F r = 0, f = 1;
    while (i)
      f /= b, r += f * (i % b), i /= b;
    return i;
  }

  constexpr F xy() { return xy(i = (i % LIM) + 1); }
};

template <typename F = float>
std::complex<F> disk(std::complex<F> xy01) noexcept {
  return std::polar<F>(std::sqrt(xy01.real()),
                       F(2) * std::numbers::pi_v<F> * xy01.imag());
}

} // namespace dyn

#endif // GRASS_HALTON_H
