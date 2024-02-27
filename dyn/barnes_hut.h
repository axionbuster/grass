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

template <class, std::unsigned_integral M = uint64_t> class View;
void run(auto const &group, auto &&process);

namespace detail {

/// A consecutive group of particles.
/// @tparam E Extra data (type of field `data`). An E object may be constructed
/// by passing two iterators to particles ("begin" and "end") or by default
/// construction, and it can be updated using the `+=` operator with an E object
template <typename E, typename I> class Group {
  template <class, std::unsigned_integral> friend class View;
  template <class A, class B> friend void run(A, B);

  /// @brief Non-empty range of particles (last is the past-the-end iterator).
  I first, last;

  /// @brief Left-child, right-sibling pointers (empty by default).
  std::shared_ptr<Group> left, right;

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

  /// Test whether this group holds exactly one particle.
  [[nodiscard]] bool single() const {
    if (first == last)
      return false;
    auto a{first};
    return ++a == last;
  }
};

} // namespace detail

/// A view to a set of particles.
/// @tparam S A set of particles (begin() and end() iterates over the particles)
/// @tparam M An unsigned integer type used for masking Morton (Z) codes.
template <class S, std::unsigned_integral M> class View {
  /// Mask (begin with the finest detail, first).
  M mask = ~M{};

  /// Iterate over particles (`begin` and `end` calls).
  S s;

  /// @brief A A group type (specialized for the type variable S).
  /// @tparam E Extra data type.
  template <class E> using GroupType = detail::Group<E, decltype(s.begin())>;

  /// Compute the groups at the current level of detail.
  /// @tparam E Extra data type.
  /// @param z Take the particle and then compute the Morton (Z) code masked by
  /// the given M-type mask.
  /// @param prior The result of this function call for one finer level of
  /// detail.
  template <class E>
  [[nodiscard]] std::shared_ptr<GroupType<E>>
  layer(auto &&z, std::shared_ptr<GroupType<E>> const &prior = {}) const {
    if (!mask || s.begin() == s.end())
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
      auto g = prior, h = prior->right;
      // Merge groups having the same Morton (Z) prefixes a and b into group q.
      // p is the root (q gets updated).
      auto p = std::make_shared<>(*g), q = p;
      // Let g be the child (left) of q.
      q->left = g;
      auto a = prefix(*(*g)->first);
      while (h) {
        auto b = prefix(*(*h)->first);
        if (a == b) {
          // If same prefix, update q:
          //  - Set boundary of past-the-last particle iterator (last).
          //  - Merge physical summary data (+=).
          q->last = h->last;
          q->data += h->data;
        } else {
          // New prefix:
          //  - New group (r).
          //  - Let r be the sibling (right) of q.
          //  - Replace the groups.
          auto r = std::make_shared<>(*h);
          q->right = r;
          std::swap(q, r);
          a = b;
          g = h;
        }
        h = h->right;
      }
      return p;
    } else {
      // First call? No problem. Turn each particle into a group.
      auto j = s.begin(), i = j++;
      auto g = std::make_shared<GroupType<E>>(i, j);
      if (i == s.end())
        return g;
      while (++i, ++j != s.end()) {
        auto h = std::make_shared<GroupType<E>>(i, j);
        h->right = g;
        std::swap(g, h);
      }
      return g;
    }
  }

public:
  View(S s) : s(s) {}

  /// Make the level of detail one level coarser.
  void coarser() {
    // Let `mask` run off to zero if at last level of detail.
    mask <<= 2;
  }

  template <class E>
  [[nodiscard]] std::shared_ptr<GroupType<E>> tree(auto &&z) const {
    auto siblings = [](auto &&l) {
      size_t c{};
      while (l)
        ++c, l = l->right;
      return c;
    };
    auto l = layer(z), c = siblings(l);
    while (l) {
      coarser(), l = layer(z, l);
      if (auto d = siblings(l); c == d)
        return l;
      else
        c = d;
    }
    return l;
  }
};

void run(auto const &group, auto &&process) {
  if (!group)
    return;
  // Apply dfs with the stack (s) of groups.
  auto s = std::vector{group};
  while (!s.empty()) {
    if (auto top = s.back(); s.pop(), !process(top) && top->left)
      // False-y value from `process`: Further detail required.
      // Push child of top.
      s.push_back(top->left);
    // Push siblings of top.
    for (auto g = group->right; !g.empty(); g = g->right)
      s.push_back(g);
  }
}

} // namespace dyn::bh32

#endif // GRASS_BARNES_HUT_H
