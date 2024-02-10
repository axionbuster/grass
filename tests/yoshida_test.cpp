#include "gtest/gtest.h"

#include <complex>

#include <yoshida.h>

using namespace std::literals::complex_literals;

TEST(YoshidaSuite, Circle0) {
  dyn::Yoshida<> yoshi{1.0f, 1.0if};
  auto accel = [](auto xy) {
    auto r = 1 / std::abs(xy);
    return -r * r * r * xy;
  };
  auto constexpr dt = 0.03125f;
  auto constexpr STEPS = 2'500'000;
  for (long i = 0; i < STEPS; i++)
    yoshi.step(dt, accel);
  auto r = std::abs(yoshi.y0),
       v = std::abs(yoshi.y1);
  auto y0 = yoshi.y0, y1 = yoshi.y1;
  auto dot = y0.real()*y1.real()+y0.imag()*y1.imag();
  ASSERT_NEAR(1.0f, r, 0.01f);
  ASSERT_NEAR(1.0f, v, 0.01f);
  ASSERT_NEAR(0.0f, dot, 0.01f);
}
