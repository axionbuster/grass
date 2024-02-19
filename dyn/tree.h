#ifndef GRASS_TREE_H
#define GRASS_TREE_H

#include <algorithm>
#include <complex>
#include <ranges>
#include <vector>

namespace dyn::tree32 {

/// @brief A node in the Barnes-Hut tree at a given depth (span zero or more
/// particles).
template <std::ranges::range R> struct Node {
  R range{};
  std::complex<float> center, centroid;
  float radius{};
};

/// @brief Given a Z-sorted (Morton-ordered) range, and a function to get the
/// prefix (at an appropriate level of detail) of each particle, as well as
/// functions to get the position and mass of said particle, call the
/// `with_node` callback whenever a new node is constructed.
template <std::ranges::range R, std::ranges::range S>
void group(R range, auto const &get_z_masked, auto const &get_xy,
           auto const &get_mass, auto const &make_range,
           auto const &with_node) {
  auto head_z = get_z_masked(*range.begin());
  for (;;) {
    auto [first, last] = std::ranges::equal_range(range.begin(), range.end(),
                                                  head_z, {}, get_z_masked);
    if (first == last)
      break;
    auto center = std::complex{0.0f, 0.0f}, centroid = center;
    auto radius = 0.0f, mass = 0.0f;
    auto count = 1.0f;
    for (auto i = first; i != last; i++) {
      auto xy = get_xy(*i);
      auto m = get_mass(*i);
      center += (xy - center) / count++;
      centroid += xy * m;
      mass += m;
    }
    centroid /= mass;
    for (auto i = first; i != last; i++) {
      auto xy = get_xy(*i);
      radius = std::max(radius, xy - center);
    }
    with_node(Node<S>{make_range(first, last), center, centroid, radius});
    if (last == range.end())
      break;
    head_z = get_z_masked(*last);
  }
}

} // namespace dyn::tree32

#endif //GRASS_TREE_H
