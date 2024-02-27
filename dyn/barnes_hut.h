#ifndef GRASS_BARNES_HUT_H
#define GRASS_BARNES_HUT_H

#include <algorithm>
#include <array>
#include <complex>
#include <concepts>
#include <cstdint>
#include <memory>
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

namespace detail {

template <class, std::unsigned_integral M = uint64_t> class View;
void run(auto const &group, auto &&process);
auto tree(auto &&, auto &&, auto &&);

/// A consecutive group of particles.
/// @tparam E Extra data (type of field `data`). An E object may be constructed
/// by passing two iterators to particles ("begin" and "end") or by default
/// construction, and it can be updated using the `+=` operator with an E object
template <typename E, typename I> class Group {
  template <class, std::unsigned_integral> friend class View;
  friend void run(auto const &, auto &&);
  friend auto tree(auto &&, auto &&, auto &&);

  /// @brief Non-empty range of particles (last is the past-the-end iterator).
  I first, last;

  /// @brief Left-child, right-sibling pointers (empty by default).
  std::shared_ptr<Group> child, sibling;

  /// @brief Extra data.
  E data{};

  /// Create a group of particles using iterators to particles.
  Group(I first, I last) : first{first}, last{last}, data{first, last} {}

  /// Iterate over the particles.
  [[nodiscard]] I begin() { return first; }

  /// Iterate over the particles.
  [[nodiscard]] I end() { return last; }

  /// Iterate over the particles.
  [[nodiscard]] I begin() const { return first; }

  /// Iterate over the particles.
  [[nodiscard]] I end() const { return last; }

public:
  /// Test whether this group holds exactly one particle.
  [[nodiscard]] bool single() const {
    if (first == last)
      return false;
    auto a{first};
    return ++a == last;
  }

  /// Return a reference to the data.
  [[nodiscard]] E &operator()() { return data; }

  [[nodiscard]] E const &operator()() const { return data; }
};

/// A view to a set of particles.
/// @tparam I Iterator of particles.
/// @tparam M An unsigned integer type used for masking Morton (Z) codes.
template <class I, std::unsigned_integral M> class View {
  friend void run(auto const &, auto &&);
  friend auto tree(auto &&, auto &&, auto &&);

  /// Mask (begin with the finest detail, that is, all ones, first).
  /// At all zeroes, stop processing.
  M mask = ~M{};

  /// Iterators over the particles. Possibly empty.
  I first, last;

  /// Compute the groups at the current level of detail.
  /// @tparam E Extra data type.
  /// @param z Take the particle and then compute the Morton (Z) code masked by
  /// the given M-type mask.
  /// @param prior The result of this function call for one finer level of
  /// detail (see `coarser`).
  template <class E>
  [[nodiscard]] std::shared_ptr<Group<E, I>>
  layer(auto &&z, std::shared_ptr<Group<E, I>> const &prior = {}) const {
    if (!mask || first == last)
      // Empty output is used as a sentinel to halt processing by the free
      // function `run`.
      return {};

    // prefix: Compute the masked Morton (Z) code, which computes the prefix
    // bits (less significant bits cut off).
    std::unsigned_integral auto m = mask;
    auto prefix = [&z, m](auto &&p) { return z(p, m); };

    if (prior) {
      // One finer level of detail exists in `prior`.
      // Merge the groups, then, instead of recalculating everything.
      // Two-pointer solution (let g and h be pointers to groups).
      auto g = prior;
      // Merge groups having the same Morton (Z) prefixes a and b into group q.
      // p is the root (q gets updated).
      auto p = std::make_shared<>(*g), q = p;
      // Let g be the child of q.
      q->child = g;
      // Let a and b be the Morton (Z) prefixes of g and h.
      auto a = prefix(*(*g)->first); // first particle.
      while ((g = g->sibling)) {
        auto b = prefix(*(*g)->first); // (ditto).
        if (a == b) {
          // If same prefix, update q:
          //  - Set boundary of past-the-last particle iterator (last).
          //  - Merge physical summary data (+=).
          q->last = g->last;
          q->data += g->data;
        } else {
          // New prefix:
          //  - New group (r).
          //  - Let r remember where it comes from (g).
          //  - Let r be the sibling of q.
          //  - Replace the groups.
          auto r = std::make_shared<>(*g);
          r->child = g;
          q->sibling = r;
          std::swap(q, r);
          a = b;
        }
      }
      return p;
    } else {
      // First call? No problem. Turn each particle into a group.
      // Let i and j be iterators to the particles.
      // Know that first != last here.
      auto j = first, i = j++;
      auto g = std::make_shared<Group<E, I>>(i, j);
      if (j == last)
        return g;
      while (++i, ++j != last) {
        auto h = std::make_shared<Group<E, I>>(i, j);
        h->sibling = g;
        std::swap(g, h);
      }
      return g;
    }
  }

  View(I first, I last) : first{first}, last{last} {}

  /// Make the level of detail one level coarser.
  void coarser() {
    // Let `mask` run off to zero if at last level of detail.
    mask <<= 2;
  }

  /// Construct a tree and then return a root node.
  template <class E>
  [[nodiscard]] std::shared_ptr<Group<E, I>> tree(auto &&z) const {
    auto siblings = [](auto &&l) {
      size_t c{};
      while (l)
        ++c, l = l->sibling;
      return c;
    };
    auto l = layer(z), c = siblings(l);
    while (l) {
      coarser(), l = layer(z, l);
      // Count of groups decreases if and only if further processed
      // (otherwise idempotent).
      if (auto d = siblings(l); c == d)
        // Cut tree to prevent duplicated visits during depth-first search.
        return l;
      else
        c = d;
    }
    return l;
  }
};

/// Construct a tree given an iterator to particles.
/// @param first Beginning iterator to the particles.
/// @param last Ending iterator to the particles (possibly equal to first).
/// @param z Given a particle and an unsigned integer mask, compute the Morton
/// code (Z code) as an unsigned integer and then apply the mask by bitwise AND.
auto tree(auto &&first, auto &&last, auto &&z) {
  auto v = View{first, last};
  return v.tree(z);
}

/// Run the group.
/// @param group Return value of a call to `tree`.
/// @param process Given a reference to a group,
void run(auto const &group, auto &&process) {
  if (!group)
    return;
  // Apply depth-first search with the stack (s) of groups.
  auto s = std::vector{group};
  while (!s.empty()) {
    if (auto top = std::move(s.back()); s.pop(), !process(top) && top->child)
      // False-y value from `process`: Further detail required.
      // Push child of top.
      s.push_back(top->child);
    // Push siblings of top.
    for (auto g = group->sibling; !g.empty(); g = g->sibling)
      s.push_back(g);
  }
}

} // namespace detail

using detail::run;
using detail::tree;

} // namespace dyn::bh32

#endif // GRASS_BARNES_HUT_H
