#ifndef GRASS_BARNES_HUT_H
#define GRASS_BARNES_HUT_H

#include <algorithm>
#include <array>
#include <complex>
#include <concepts>
#include <cstdint>
#include <optional>
#include <vector>

namespace dyn::bh32 {

namespace detail {

/// @brief Interleave the bits of two 32-bit words (re, im) so that the word
/// im is placed at the odd-numbered bits (including the most significant bit)
/// while re is placed at the even-numbered bits (including the least
/// significant bit). Bits are numbered from 0 (LSB) to 61 (MSB) inclusive.
constexpr uint64_t interleave32(uint32_t re, uint32_t im) {
  // Modification of "Interleave by Binary Magic Numbers"
  // (http://graphics.stanford.edu/~seander/bithacks.html#InterleaveBMN)
  struct Help {
    uint64_t mask, shift;
  };

  // 32-bit words. Hence, 5.
  // 0101..., 00110011...,
  // 00001111..., 0000000011111111...,
  // including until the entire word is like 0...0_1...1_0...0_1...1.
  // (Reverse order, however.)
  std::array<Help, 5> constexpr H{
      Help{0x0000ffff0000ffff, 16}, Help{0x00ff00ff00ff00ff, 8},
      Help{0x0f0f0f0f0f0f0f0f, 4}, Help{0x3333333333333333, 2},
      Help{0x5555555555555555, 1}};

  // Zero-extend each word.
  std::array<uint64_t, 2> W{re, im};

  // Spread the bits out with zeros in between.
  // (ex) w = 0b1011 -> 0b01'00'01'01
  for (auto &&w : W)
    for (auto &&h : H)
      w = (w | (w << h.shift)) & h.mask;

  // Imaginary first.
  return W[0] | (W[1] << 1);
}

/// @brief Group consecutive particles with the same Z-code given a Z-sorted
/// range of particles. (A Z-code is also known as a Morton code).
/// @param begin Starting iterator for the particles.
/// @param end Past-the-end iterator for the particles, possibly equal to begin.
/// @param z Given what is like a const lvalue reference to a particle, compute
/// its z-code, perhaps masked (low bits flushed to zero). The return values of
/// `z` satisfy this contract: If `a` is returned by a call to `z`, then a == a.
/// @param g Callback with two arguments, both iterators to the particles, to
/// be called when a range of particles is announced.
void group(auto begin, auto const end, auto &&z, auto &&g) {
  while (begin != end) {
    // Apply binary search for the consecutive range.
    // (Same Z-code as *begin).
    auto [first, last] = std::ranges::equal_range(begin, end, z(*begin), {}, z);

    // Assume that first != last due to the contract fulfilled by `equal_range`,
    // and the reflexivity of equality comparison (==) of the return values of
    // `z` (an important precondition).

    // Reset `begin` and then announce it (g).
    g(first, begin = last);
  }
}
} // namespace detail

/// @brief Compute the Morton (Z) code of a complex number xy assuming a squared
/// grid by pre-multiplying the factor `Precision` to each component of xy.
/// If either component is not a finite number after the scaling, the result
/// is {} (no value). It is recommended to cache the answer because of the
/// typically high overhead of the computation.
template <uint32_t Precision = 512>
std::optional<uint64_t> morton(std::complex<float> xy) {
  xy *= float(Precision);
  // Use a strict inequality because float(INT32_MAX) is actually greater than
  // INT32_MAX (when both are cast to double).
  if (std::abs(xy.real()) < float(INT32_MAX) &&
      std::abs(xy.imag()) < float(INT32_MAX)) {
    // Flip sign bit before bit-casting to unsigned to preserve order.
    auto sgn = uint32_t(0x8000'0000ul);
    auto x = int32_t(xy.real()) ^ sgn; // promoted to unsigned.
    auto y = int32_t(xy.imag()) ^ sgn; // (ditto).
    return detail::interleave32(x, y);
  }
  return {};
}

/// A consecutive group of particles.
/// @tparam E Extra data (type of field `data`). An E object may be constructed
/// by passing two iterators to particles ("begin" and "end") or by default
/// construction, and it can be updated using the `+=` operator with an E object
template <typename E, typename I> struct Group {
  I first, last;
  E data{};
  /// Create a group of particles using iterators to particles.
  Group(I first, I last) : first{first}, last{last}, data{first, last} {}
  /// Summarize a range of smaller groups.
  Group(auto g_first, auto g_last) {
    first = g_first->first;
    while (g_first != g_last)
      last = g_first->last, data += g_first++->data;
  }
  /// Iterate over the particles.
  [[nodiscard]] I begin() { return first; }
  /// Iterate over the particles.
  [[nodiscard]] I end() { return last; }
  /// Iterate over the particles.
  [[nodiscard]] I begin() const { return first; }
  /// Iterate over the particles.
  [[nodiscard]] I end() const { return last; }
  /// Test whether this group holds exactly one particle.
  [[nodiscard]] bool single() const {
    if (first == last)
      return false;
    auto a{first};
    return ++a == last;
  }
  /// Test: Are the iterators equal (disregarding `data`)?
  [[nodiscard]] bool operator==(Group const &g) const {
    // Ignore extra data.
    return this == &g || (first == g.first && last == g.last);
  }
};

/// A view to a set of particles.
/// @tparam S A set of particles (begin() and end() iterates over the particles)
/// @tparam M An unsigned integer type used for masking Morton (Z) codes.
template <class S, std::unsigned_integral M = uint64_t> class View {
  /// Mask (begin with the finest detail, first).
  M mask = ~M{};

  /// Iterate over particles (`begin` and `end` calls).
  S s;

public:
  View(S s) : s(s) {}

  /// Return type of `groups`.
  template <class E> using Groups = std::vector<Group<E, decltype(s.begin())>>;

  /// Compute the groups at the current level of detail.
  /// @param z Take the particle and then compute the Morton (Z) code masked by
  /// the given M-type mask.
  /// @param prior The result of this function call for one finer level of
  /// detail.
  template <class E>
  [[nodiscard]] Groups<E> groups(auto &&z, Groups<E> const &prior = {}) const {
    if (!mask || s.begin() == s.end())
      // Empty output used as sentinel to halt processing by the free functions:
      //  - `levels`
      //  - `run`
      return {};

    // prefix: Compute the masked Morton (Z) code, which computes the prefix
    // bits (less significant bits cut off).
    auto m = mask;
    auto prefix = [&z, m](auto &&p) { return z(p, m); };

    if (prior.empty()) {
      // First call? No problem. Construct the groups.
      Groups<E> novel{};
      auto push = [&novel](auto f, auto l) { novel.emplace_back(f, l); };
      detail::group(s.begin(), s.end(), prefix, push);
      return novel;
    }

    // One finer level of detail exists in `prior`.
    Groups<E> novel{};

    // Merge the groups, then, instead of recalculating everything.
    // Two-pointer solution: g and h are iterators to the groups in `prior`.
    auto g = prior.begin(), h = g;
    // Merge groups having the same Morton (Z)-prefixes (a, b).
    auto a = prefix(*g->first);
    while (++h != prior.end()) {
      auto b = prefix(*h->first);
      // If same prefix, don't do anything specific.
      // If new prefix, create a group by merging many groups.
      if (a != b) {
        // Yes, new prefix. Treat h as past-the-end iterator to the range of
        // groups. Pass g and h to the "range of groups" constructor of Group.
        novel.emplace_back(g, h);
        g = h;
        a = b;
      }
    }
    // Handle runoff.
    if (g != h)
      novel.emplace_back(g, h);
    return novel;
  }

  /// Make the level of detail one level coarser.
  void coarser() {
    // Let `mask` run off to zero if at last level of detail.
    mask <<= 2;
  }
};

/// Compute all levels given a view and a way to compute the Morton codes given
/// a particle and a mask (see View::groups for information on `z`).
template <class E> auto levels(auto &&view, auto &&z) {
  // (Typing hacks)
  auto a = view.template groups<E>(z);
  auto l = std::vector<decltype(a)>{std::move(a)};
  do
    // Regurgitate.
    view.coarser(), l.push_back(view.groups(z, l.back()));
  while (!l.back().empty());
  l.pop_back();
  // Cull consecutive layers with the identical groups.
  l.erase(std::unique(l.begin(), l.end()), l.end());
  return l;
}

/// Process the levels.
/// @param levels What is returned by the free function `levels`.
/// @param process Return nonzero if a group should be discarded; zero if it
/// should be broken up into finer pieces.
void run(auto const &levels, auto &&process) {
  if (levels.rbegin() == levels.rend())
    return;

  // Below, `novel` and `prior` are collections of groups.
  // A "group" is a range of particles (plus any extra data).

  // Filter `novel` for inclusion in groups in `prior` knowing that the
  // groups are sorted in Morton order.
  auto begin = levels.rbegin();
  auto prior{*begin};
  while (++begin != levels.rend()) {
    auto const &novel = *begin;
    // copy: Filtered copy of `novel`.
    decltype(prior) copy;
    // Iterators to a prior group (pg) and a novel group (ng).
    auto pg = prior.begin();
    auto ng = novel.begin();
    if (ng == novel.end() || pg == prior.end())
      break;

    // Skip.
    while (ng->begin() != pg->begin())
      ++ng;

    // Cull.
    while (pg != prior.end()) {
      if (ng->begin() == pg->begin()) {
        while (ng->end() != pg->end())
          copy.push_back(*ng++);
        copy.push_back(*ng++), ++pg;
      } else
        while (ng->begin() != pg->begin())
          ++ng;
    }

    // Process the chosen groups in `copy`.
    std::erase_if(copy, process);

    // Break if no more.
    if (copy.empty())
      break;

    // Learn which branches are kept and which are pruned.
    prior = copy;
  }
}

} // namespace dyn::bh32

#endif // GRASS_BARNES_HUT_H
