#include "gtest/gtest.h"

#include <complex>

#include <circle.h>

/// @brief See `CircleTest0.png` for reference.
class CircleTest0 : public testing::Test {
protected:
  const dyn::Circle<float> circle{0, 2.4f};
  // Rectangle's less-less corner
  const std::complex<float> less_less{-4.0f, -4.0f};
  // Rectangle's greater-greater corner
  const std::complex<float> greater_greater{-2.0f, -2.0f};
};

TEST_F(CircleTest0, Disjoint0) {
  ASSERT_FALSE(dyn::disk_rect_intersect(circle, less_less, greater_greater));
}
