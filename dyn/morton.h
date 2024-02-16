#ifndef GRASS_MORTON_H
#define GRASS_MORTON_H

#include <array>
#include <bit>
#include <complex>
#include <cstdint>

namespace dyn {

namespace morton::detail {
constexpr uint32_t order32(float x) {
  // Due to Tropf (2021).
  auto sgn = uint32_t(1) << 31, i = std::bit_cast<uint32_t>(x);
  return (i & sgn) ? ~i : (i | sgn);
}

constexpr uint64_t interleave32(uint32_t re, uint32_t im) {
  // Modification of "Interleave by Binary Magic Numbers"
  // (http://graphics.stanford.edu/~seander/bithacks.html#InterleaveBMN)
  struct Help {
    // 0101..., 00110011...,
    // 00001111..., 0000000011111111...,
    // including until the entire word is like 0...0_1...1_0...0_1...1
    // (In `H` below, the above is done in reverse order.)
    uint64_t mask;
    // Number of bits to shift.
    uint64_t shift;
  };

  // 32-bit words. Hence, 5.
  std::array<Help, 5> constexpr H{
      Help{0x0000ffff0000ffff, 16}, Help{0x00ff00ff00ff00ff, 8},
      Help{0x0f0f0f0f0f0f0f0f, 4}, Help{0x3333333333333333, 2},
      Help{0x5555555555555555, 1}};

  std::array<uint64_t, 2> W{re, im};
  for (auto &&w : W)
    for (auto &&h : H) {
      // "Duplicate" the bits (shift-then-OR) and then mask result.
      w = (w | (w << h.shift)) & h.mask;
    }
  // Imaginary first.
  return W[0] | (W[1] << 1);
}
} // namespace morton::detail

constexpr uint64_t morton32(std::complex<float> xy) {
  using namespace morton::detail;
  return interleave32(order32(xy.real()), order32(xy.imag()));
}

} // namespace dyn

#endif // GRASS_MORTON_H
