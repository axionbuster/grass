#ifndef GRASS_BARNES_HUT_H
#define GRASS_BARNES_HUT_H

#include <algorithm>
#include <complex>
#include <concepts>
#include <ranges>
#include <vector>

namespace dyn::bh32 {

/// @brief Given a Z-sorted (Morton-ordered) range and a function to get the
/// prefix (at an appropriate level of detail) of the Morton code of each
/// particle, construct an array of Node objects at that level of detail.
void group(auto begin, auto const end, auto &&z, auto &&node) {
  while (begin != end) {
    auto [first, last] = std::ranges::equal_range(begin, end, z(*begin), {}, z);
    if (first != last)
      node(first, last), begin = last;
    else
      ++begin;
  }
}

void dfs(auto begin, auto const end, auto &&far, auto &&react) {
  if (begin != end)
    do
      if (far(*begin))
        react(*begin);
      else
        dfs(begin->kids.begin(), begin->kids.end(), std::forward(far),
            std::forward(react));
    while (++begin != end);
}

} // namespace dyn::bh32

#endif // GRASS_BARNES_HUT_H
