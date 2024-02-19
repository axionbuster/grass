#include "gtest/gtest.h"

#include <algorithm>
#include <morton.h>
#include <random>
#include <tree.h>

TEST(Trees, Trees0) {
  // Agree with Morton order?
  std::mt19937 rng(1234);
  std::uniform_real_distribution<float> dist(-2.0f, 2.0f);
  std::vector<std::complex<float>> pp;
  pp.reserve(32767);
  for (int i = 0; i < 32767; i++)
    pp.emplace_back(dist(rng), dist(rng));
  std::ranges::sort(pp.begin(), pp.end(), {}, dyn::fixedmorton32<512>);
  auto constexpr PRECISION = 3;
  auto constexpr LL = std::complex{-1.0f, -1.0f};
  auto constexpr SIDE = 2.0f;
  auto get_xy = [&pp](size_t i) { return pp[i]; };
  auto get_mass = [](size_t _) { return 1.0f; };
  auto tree = dyn::tree32::Tree::from_z_ordered(PRECISION, pp.size(), LL, SIDE,
                                                get_xy, get_mass);
  // FIXME: Currently a smoke test, but I will add an invariant check soon.
}
