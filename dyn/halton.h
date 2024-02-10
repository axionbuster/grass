#ifndef GRASS_HALTON_H
#define GRASS_HALTON_H

#include <complex>
#include <numbers>

namespace dyn {

template <typename F = float, unsigned short b = 2, unsigned int LIM = 0x1000>
class Halton {
  unsigned short i{};

public:
  constexpr static F x01(unsigned short i) {
    F r = 0, f = 1;
    while (i)
      f /= b, r += f * (i % b), i /= b;
    return i;
  }

  constexpr F x01() { return x01(i = (i % LIM) + 1); }
};

} // namespace dyn

#endif // GRASS_HALTON_H
