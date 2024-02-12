#ifndef GRASS_IRHALL_H
#define GRASS_IRHALL_H

// Irwin-Hall approximation to the uniform distribution.

namespace phy {

/// @brief Approximate a normal distribution using a function that generates a
/// uniform (0,1) variate each time it is called using the Irwin-Hall
/// approximation (and the central limit theorem).
template <typename F = float, typename A> F irhnormal(A u01) {
  F s{};
  for (unsigned char i = 0; i < 12; i++)
    s += u01();
  return s - F(6);
}

} // namespace phy

#endif // GRASS_IRHALL_H
