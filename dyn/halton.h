#ifndef GRASS_HALTON_H
#define GRASS_HALTON_H

/// @file halton.h
/// @brief The Halton low-discrepancy sequence.

namespace dyn {

/// @brief The Halton low-discrepancy sequence. A low-discrepancy sequence
/// generates points on the unit interval (0, 1) as evenly as possible for any
/// number of samples.
/// @tparam F A floating-point type.
/// @tparam b "Base" (a small prime number).
/// @tparam LIM Largest possible value of "index" (inclusive).
template <typename F = float, short b = 2, unsigned int LIM = 0x1000>
class Halton {
  /// @brief Last index.
  short i{};

public:
  /// @brief Generate a number in the open interval (0, 1).
  /// @param i Index, not zero.
  constexpr static F x01(short i) {
    auto r = F(0), f = F(1);
    while (i)
      f /= b, r += f * (i % b), i /= b;
    return r;
  }

  /// @brief Generate a number in the open interval (0, 1) and advance the
  /// internal index, wrapping inclusively around the specified limit (see `LIM`
  /// template parameter).
  constexpr F x01() { return x01(i = (i % LIM) + 1); }
};

} // namespace dyn

#endif // GRASS_HALTON_H
