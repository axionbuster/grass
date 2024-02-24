#include "gtest/gtest.h"

#include <complex>

#include <circle.h>

/// @brief See `CircleTest0.png` for reference.
class CircleTest0 : public testing::Test {
protected:
  // xy, r; rectangles (2 corners)
  dyn::Circle<float> const circle{0, 2.4f};
  std::complex<float> const less_less{-4.0f, -4.0f};
  std::complex<float> const greater_greater{-2.0f, -2.0f};
};

TEST_F(CircleTest0, Disjoint0) {
  ASSERT_FALSE(dyn::disk_arrect_isct(circle, less_less, greater_greater));
}

class CircleTest1 : public testing::Test {
protected:
  dyn::Circle<float> const circle{0, 2.0f};
  std::complex<float> const less_less{-5.0f, 1.0f},
      greater_greater{-1.0f, 5.0f};
};

TEST_F(CircleTest1, In0) {
  ASSERT_TRUE(dyn::disk_arrect_isct(circle, less_less, greater_greater));
}

class CircleTest2 : public testing::Test {
protected:
  dyn::Circle<float> const circle{0, 3.0f};
  std::complex<float> const less_less{0.6f, 2.7f}, greater_greater{1.8f, 3.6f};
};

TEST_F(CircleTest2, In0) {
  ASSERT_TRUE(dyn::disk_arrect_isct(circle, less_less, greater_greater));
}

class CircleTest3 : public testing::Test {
protected:
  dyn::Circle<float> const circle{0, 4.0f};
  std::complex<float> const c{-2.0f, -2.0f}, e{0, -1.0f}, l{-2.0f, -5.0f},
      n{5.0f, 3.0f}, s{7.0f, 2.0f}, u{9.0f, 3.0f};
};

TEST_F(CircleTest3, In0) { ASSERT_TRUE(dyn::disk_arrect_isct(circle, c, e)); }

TEST_F(CircleTest3, In1) { ASSERT_TRUE(dyn::disk_arrect_isct(circle, l, n)); }

TEST_F(CircleTest3, Out0) { ASSERT_FALSE(dyn::disk_arrect_isct(circle, s, u)); }
