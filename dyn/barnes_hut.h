#ifndef GRASS_BARNES_HUT_H
#define GRASS_BARNES_HUT_H

#include <algorithm>

namespace dyn::bh32 {

/// @brief Given a Z-sorted (Morton-ordered) range and a function to get the
/// prefix (at an appropriate level of detail) of the Morton code of each
/// particle, construct an array of Node objects at that level of detail.
void group(auto begin, auto const end, auto &&z, auto &&grp) {
  while (begin != end) {
    auto [first, last] = std::ranges::equal_range(begin, end, z(*begin), {}, z);
    if (first != last)
      grp(first, begin = last);
    else
      ++begin;
  }
}

} // namespace dyn::bh32

#endif // GRASS_BARNES_HUT_H
