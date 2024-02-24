#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <barnes_hut.h>
#include <cstdint>
#include <iomanip>
#include <iostream>

TEST(MortonDetail, Interleave0) {
  auto a = uint32_t(0xffffffff), b = uint32_t(0x00000000);
  auto e = uint64_t(0x55555555'55555555),
       g = dyn::bh32::detail::interleave32(a, b);
  ASSERT_EQ(e, g) << "hex (e) = 0x" << std::hex << std::setw(16) << e
                  << "\nhex (g) = 0x" << std::setw(16) << g;
}

TEST(Morton, Fixed512_0) {
  auto X = 4194304; // INT32_MAX / 512
  std::complex<float> a{float(X), float(X)};
  auto b = dyn::bh32::morton<512>(a);
  ASSERT_FALSE(b.has_value());
}

TEST(Morton, Fixed512_1) {
  auto X = 12345;
  std::complex<float> a{float(X), float(X)};
  auto b = dyn::bh32::morton<512>(a);
  ASSERT_TRUE(b.has_value());
}

TEST(Morton, Fixed512_2) {
  typedef std::complex<float> C;
  // Already sorted in Z-order (Morton order; no change expected).
  std::array<C, 4> z_in{C{-12.0f, -11.0f}, C{24.0f, -3.23f}, C{-11.0f, 4.8f},
                        C{1.2f, 3.4f}};
  auto z_out{z_in};
  std::ranges::sort(z_out.begin(), z_out.end(), {}, dyn::bh32::morton<512>);
  ASSERT_EQ(z_in, z_out);
}

TEST(Morton, Fixed512_3) {
  typedef std::complex<float> C;
  // In incorrect Z-order (must sort).
  std::array<C, 2> z_in{C{11.0f, 3.3f}, {-2.0f, 0.2f}};
  std::array<C, 2> z_expect{C{-2.0f, 0.2f}, {11.0f, 3.3f}};
  auto z_out{z_in};
  std::ranges::sort(z_out.begin(), z_out.end(), {}, dyn::bh32::morton<512>);
  ASSERT_EQ(z_expect, z_out);
}

// Notes: Order of NaN and infinite values are unspecified.
