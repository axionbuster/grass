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

// Used by `run_level`.
enum {
  /// @brief Keep the particles under this group.
  KEEP [[maybe_unused]],
  /// @brief Remove the particles under this group.
  REMOVE,
};

/// @brief Run through a "level" (a collection of groups) in reverse. Each
/// group contains a non-empty range of particles.
/// @param g_rbegin Group iterator begin (reverse).
/// @param g_rend Group iterator end (reverse).
/// @param process Process a group and decide whether to keep its particles
/// (KEEP) or remove them (REMOVE).
/// @param p_remove_range Remove a non-empty range of particles.
void run_level(auto g_rbegin, auto const g_rend, auto &&process,
               auto &&p_remove_range) {
  if (g_rbegin != g_rend)
    for (auto i = g_rbegin; i != g_rend; i++)
      if (process(*i) == REMOVE)
        p_remove_range(i->begin(), i->end());

  // Original from Notes:

  // Check for stopping conditions (no next level, no particles).
  //  "no next level" -> Done.
  //  "no particles" -> N/A.
  // Build the group buffer, applying the mask. -> N/A.
  // Process from right to left, mark for removal as needed. -> Done.
  // If marked for removal, remove all particles in the range. -> Done.
  // Continue to next group if exists; otherwise, do following: -> N/A.
  // Discard the group buffer. -> N/A.
  // Empty the particle vector. -> N/A.
  // Yield (must go to the next level). -> N/A.
}

} // namespace dyn::bh32

#endif // GRASS_BARNES_HUT_H
