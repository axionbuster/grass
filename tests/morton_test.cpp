#include "gtest/gtest.h"

#include <cstdint>
#include <iomanip>
#include <iostream>
#include <morton.h>

TEST(MortonDetail, Interleave0) {
  auto a = uint32_t(0xffffffff), b = uint32_t(0x00000000);
  auto e = uint64_t(0x55555555'55555555),
       g = dyn::morton::detail::interleave32(a, b);
  ASSERT_EQ(e, g) << "hex (e) = 0x" << std::hex << std::setw(16) << e
                  << "\nhex (g) = 0x" << std::setw(16) << g;
}
