#ifndef GRASS_BARNES_HUT_H
#define GRASS_BARNES_HUT_H

#include <array>
#include <cassert>
#include <complex>
#include <cstdint>
#include <memory>
#include <optional>
#include <stack>
#include <utility>

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

template <class, class I> auto tree(I, I, auto &&);
void run(auto const *, auto &&);
struct DeleteGroup;

template <class I, class E> class Group {
  template <class, class J> friend auto tree(J, J, auto &&);
  friend void run(auto const *, auto &&);
  friend struct DeleteGroup;

  I first, last;
  Group *child{}, *sibling{};

public:
  E extra;

  ~Group() = default;

private:
  Group(Group const &g) : first{g.first}, last{g.last}, extra{g.extra} {}
  Group(I first, I last) : first{first}, last{last}, extra{first, last} {}
  Group(Group *a, Group *b) : first{a->first}, last{b->last} {
    child = a;
    extra = a->extra;
    for (Group *c = a->sibling; c && c != b; c = c->sibling)
      extra += c->extra;
  }
  Group(Group &&g) noexcept
      : first{std::move(g.first)}, last{std::move(g.last)},
        extra{std::move(g.extra)} {}
  Group &operator=(Group const &g) {
    if (this != &g)
      first = g.first, last = g.last, extra = g.extra;
    return *this;
  }
  Group &operator=(Group &&g) noexcept {
    if (this != &g) {
      first = std::move(g.first), last = std::move(g.last),
      extra = std::move(g.extra);
    }
    return *this;
  }

  static void depth_first_delete(Group *g) {
    if (!g)
      return;
    std::stack<Group *> v;
    while (!v.empty()) {
      auto h = v.top();
      v.pop();
      for (auto a = h->child; a; a = a->sibling)
        v.push(a);
      delete h;
    }
  }

  static void depth_first(Group const *g, auto &&deeper) {
    if (!g)
      return;
    std::stack<Group const *> v;
    while (!v.empty()) {
      auto h = v.top();
      v.pop();
      for (auto i = h->child; i; i = i->sibling)
        if (deeper(i->extra))
          v.push(i);
    }
  }
};

struct DeleteGroup {
  template <class G> void operator()(G *g) const { G::depth_first_delete(g); }
};

inline DeleteGroup constexpr deleteGroup;

template <class E, class I> auto tree(I const first, I const last, auto &&z) {
  struct {
    uint64_t mask = ~uint64_t{};
    void shift() { mask <<= 2; }
  } state;
  using G = Group<I, E>;
  using P = std::shared_ptr<G>;

  // Turn every particle into a group.
  if (first == last)
    // Degeneracy 0 (no particles).
    return P{};
  // Let a and b be iterators to the particles exactly one particle apart.
  I b = first, a = b++;
  auto *g = new G{a, b};
  if (b == last)
    // Degeneracy 1 (one particle).
    return P{g, deleteGroup};
  // Ordinary case (two+ particles).
  // Let g and h be pointers to groups.
  auto *h = g;
  do {
    // Still, a and b are one particle apart.
    ++a, ++b;
    // g, h, and i are groups with the same level of detail.
    h->sibling = new G{a, b};
    h = h->sibling;
  } while (b != last);

  // Let's move out (lower level of detail).
  // Merge groups with the same Morton code (Z-code) prefixes.

  // h has (will have) a lower level of detail than g.
  h = new G{*g};
  h->child = g;
  state.shift();
  while (state.mask) {
    auto prefix = [&z, &state](auto &&a) { return z(a, state.mask); };
    auto c = prefix(*g->first);
    // same: Built exactly the same layer?
    auto same = true;
    for (auto *i = g->sibling; i; i = i->sibling) {
      // g (unchanging in this loop) and i have the same level of detail.
      // h has a lower level of detail than both.
      auto d = prefix(*i->first);
      if (c != d) {
        h = std::exchange(h->sibling, new G{g, i});
        g = i;
        c = d;
        same = false;
      }
    }

    // If built exactly the same layer then just use the old layer.
    g = h;
    if (same) {
      auto *j = std::exchange(g->child, {});
      deleteGroup(g);
      g = j;
    }

    // Build new layer.
    h = new G{*g};
    h->child = g;
    state.shift();
  }

  // Mark special (empty range) -> root.
  h->first = h->last;
  assert(!h->sibling);
  return P{h, deleteGroup};
}

void run(auto const *group, auto &&f) {
  std::remove_reference_t<decltype(*group)>::depth_first(group, f);
}

} // namespace detail

using detail::run;
using detail::tree;

} // namespace dyn::bh32

#endif // GRASS_BARNES_HUT_H
