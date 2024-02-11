#include "gtest/gtest.h"

#include <cmath>
#include <complex>
#include <numbers>

#include <kahan.h>
#include <newton.h>
#include <yoshida.h>

using namespace std::literals::complex_literals;

class NewtonSuite : public testing::Test {
protected:
  dyn::Circle<> const c0{0.0f + 0.0if, 0.04f};
  float const m0 = 1.0f, r1 = 0.04f;
  float const dt = 0.0625f;
  long const steps = 2'500'000;
  dyn::Gravity<float, 150> gr;

  std::complex<float> common_accel(std::complex<float> xy,
                                   dyn::Circle<float> src) {
    dyn::Circle<> c1{xy, r1};
    return gr.field(c1, src, 1.0f);
  }
};

TEST_F(NewtonSuite, YoshidaCircle0) {
  dyn::Yoshida<> yoshi{1.0f, 1.0if};

  auto accel = [&](auto xy) { return this->common_accel(xy, c0); };
  for (long i = 0; i < steps; i++)
    yoshi.step(dt, accel);

  auto r = std::abs(yoshi.y0), v = std::abs(yoshi.y1);
  ASSERT_NEAR(1.0f, r, 0.01f);
  ASSERT_NEAR(1.0f, v, 0.01f);
}

TEST_F(NewtonSuite, PassThrough0) {
  dyn::Circle<> c0{2.1f - 4.5if, 1.0f};
  auto constexpr S2 = std::numbers::sqrt2_v<float>;
  dyn::Yoshida<> yoshi{c0 + S2 + S2 * 1.0if, 0};

  // 90 steps : 1 second = 1800 steps : 20 seconds
  auto constexpr steps = 1800;
  auto accel = [&](auto xy) { return this->common_accel(xy, c0); };
  for (int i = 0; i < steps; i++) {
    yoshi.step(dt, accel);
    auto r = std::abs(yoshi.y0 - c0);
    ASSERT_LE(r, 2.05f) << "at i = " << i;
  }

  auto vsgn0 = std::signbit(yoshi.y1.real());
  for (;;) {
    yoshi.step(dt, accel);
    auto vsgn1 = std::signbit(yoshi.y1.real());
    if (vsgn0 != vsgn1)
      break;
  }

  auto r = std::abs(yoshi.y0 - c0);
  ASSERT_NEAR(2.0f, r, 0.1f);

  auto t = std::arg(yoshi.y0 - c0);
  auto constexpr mpi34 = -3 * std::numbers::pi_v<float> / 4;
  auto constexpr pi4 = std::numbers::pi_v<float> / 4;
  if (t > 0)
    ASSERT_NEAR(pi4, t, 0.05f);
  else if (t < 0)
    ASSERT_NEAR(mpi34, 4, 0.05f);
  else
    FAIL() << "(angle t = " << t << " neither positive nor negative)";
}

TEST_F(NewtonSuite, Inside0) {
  // Inside a large circle (c0), no force is felt by a test particle.

  dyn::Circle<> c0 = this->c0;
  c0.radius = 1.0f;

  // Position and velocity, resp., of c1.
  // (Since r1 is 0.04, c1 will be fully contained inside c0).
  dyn::Yoshida<> yoshi{0.25f, 0};

  // Go!
  auto accel = [&](auto xy) { return this->common_accel(xy, c0); };
  for (long i = 0; i < steps; i++) {
    yoshi.step(dt, accel);
  }

  auto v = std::abs(yoshi.y1);
  ASSERT_EQ(v, 0);
}

TEST_F(NewtonSuite, Figure8) {
  // Data from Wikipedia
  // https://en.wikipedia.org/w/index.php?title=Three-body_problem&oldid=1199934443#Special-case_solutions

  std::complex<float> _c0{-0.97000436f, 0.24308753f},
      _v0{0.4662036850f, 0.4323657300f}, _v1{-0.93240737f, -0.86473146f};

  // In the figure-8 three-body problem, the bodies don't intersect and stay far
  // away from each other. So, make the radius small enough, so they won't
  // accidentally touch each other.
  float constexpr RADIUS = 0.025f;

  // Initial conditions (copy construct them, so I can compare them).
  dyn::Yoshida<> yoshi0{_c0, _v0}, yoshi1{0, _v1}, yoshi2{-_c0, _v0};
  dyn::Yoshida<> yoshis[3] = {yoshi0, yoshi1, yoshi2};

  // Step.
  auto constexpr dt = 0.04;
  long constexpr STEPS = 158; // stop at the period t = 6.33.
  for (long s = 0; s < STEPS; s++) {
    // Particle i feels a force from the other particles j.
    for (auto i = 0; i < 3; i++) {
      auto accel = [&](auto xy) {
        dyn::Kahan<std::complex<float>> a;
        for (auto j = 0; j < 3; j++)
          if (i != j) {
            dyn::Circle<> ci{xy, RADIUS}, cj{yoshis[j].y0, RADIUS};
            // Unit mass
            a += gr.field(ci, cj, 1.0f);
          }
        return a();
      };
      yoshis[i].step(dt, accel);
    }
  }

  // They must be near each other to about one or two decimal precision
  // after the period.

  float qs[3] = {std::abs(yoshi0.y0 - yoshis[0].y0),
                 std::abs(yoshi1.y0 - yoshis[1].y0),
                 std::abs(yoshi2.y0 - yoshis[2].y0)};

  for (int i = 0; i < 3; i++)
    ASSERT_NEAR(1.0f, qs[i], 0.1f) << "(i = " << i << ")";
}
