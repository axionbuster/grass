#ifndef GRASS_TREE_H
#define GRASS_TREE_H

#include <algorithm>
#include <complex>
#include <ranges>
#include <vector>

namespace dyn::tree32 {

/// @brief A node in the Barnes-Hut tree at a given depth (span zero or more
/// particles).
template <typename I> struct Node {
  I first, last;
  std::complex<float> center, centroid;
  float radius{};
  I begin() { return first; }
  I begin() const { return first; }
  I end() { return last; }
  I end() const { return last; }
};

/// @brief Given a Z-sorted (Morton-ordered) range, and a function to get the
/// prefix (at an appropriate level of detail) of each particle, as well as
/// functions to get the position and mass of said particle, call the
/// `with_node` callback whenever a new node is constructed.
template <typename I>
void group(I begin, I const end, auto &&get_z_masked, auto &&with_node) {
  while (begin != end) {
    auto [first, last] = std::ranges::equal_range(
        begin, end, get_z_masked(*begin), {}, get_z_masked);
    if (first != last) {
      with_node(Node{first, last, {}, {}, {}});
      begin = last;
    } else {
      begin++;
    }
  }
}

} // namespace dyn::tree32

#endif //GRASS_TREE_H
