#ifndef GRASS_BARNES_HUT_H
#define GRASS_BARNES_HUT_H

#include <algorithm>
#include <complex>
#include <concepts>
#include <ranges>
#include <vector>

namespace dyn::bh32 {

/// @brief A node in the Barnes-Hut tree at a given depth (span zero or more
/// particles).
/// @tparam E Any extra data (default-constructible).
/// @tparam I Iterator type.
template <std::default_initializable E, typename I> struct Node {
  I first, last;
  E extra;
  I begin() { return first; }
  I begin() const { return first; }
  I end() { return last; }
  I end() const { return last; }
};

/// @brief Given a Z-sorted (Morton-ordered) range and a function to get the
/// prefix (at an appropriate level of detail) of the Morton code of each
/// particle, construct an array of Node objects at that level of detail.
/// @tparam E Extra data (correspond to the same parameter in struct Node).
/// @tparam I Iterator type (correspond to the same in struct Node).
template <typename E, typename I>
void group(I begin, I const end, auto &&z, auto &&node) {
  while (begin != end) {
    auto [first, last] = std::ranges::equal_range(begin, end, z(*begin), {}, z);
    if (first != last)
      node(Node<E, I>{first, last, {}}), begin = last;
    else
      begin++;
  }
}

void dfs(auto begin, auto const end, auto &&far, auto &&react) {
  if (begin != end)
    do
      if (far(*begin))
        react(*begin);
      else
        dfs(begin->node_begin(), begin->node_end(), std::forward(far),
            std::forward(react));
    while (++begin != end);
}

} // namespace dyn::bh32

#endif // GRASS_BARNES_HUT_H
