#ifndef GRASS_BARNES_HUT_H
#define GRASS_BARNES_HUT_H

#include <array>
#include <cassert>
#include <complex>
#include <cstdint>
#include <memory>
#include <optional>
#include <stack>

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

template <class, class I> auto tree(I, I, auto &&) noexcept;
void run(auto const *, auto &&);
struct DeleteGroup;

template <class I, class E> class Group {
  template <class, class J> friend auto tree(J, J, auto &&) noexcept;
  friend void run(auto const *, auto &&);
  friend struct DeleteGroup;

  I first, last;
  Group *child{}, *sibling{};

public:
  E extra;

  ~Group() = default;

private:
  Group(I const first, I const last)
      : first{first}, last{last}, extra{first, last} {}
  Group(Group *const a, Group *const b) : first{a->first} {
    // b can be nullptr.
    extra = a->extra;
    for (auto *c = a->sibling; c && c != b; c = c->sibling)
      extra += c->extra, last = c->last;
  }
  void clear_relationships() { child = {}, sibling = {}; };

  static void depth_first_delete(Group *const g) {
    if (!g)
      return;
    size_t debug_count{};
    std::stack<Group *> v;
    v.push(g);
    while (!v.empty()) {
      debug_count++;
      auto h = v.top();
      assert(h);
      v.pop();
      for (auto a = h->child; a; a = a->sibling)
        v.push(a);
      delete h;
    }
    (void)debug_count;
  }

  /// Apply depth-first traversal. If `deeper` suggests going deeper (true),
  /// go deeper.
  static void depth_first(Group const *g, auto &&deeper) {
    if (!g)
      return;
    size_t debug_count{};
    std::stack<Group const *> v;
    v.push(g);
    while (!v.empty()) {
      auto h = v.top();
      v.pop();
      for (auto i = h->child; i; i = i->sibling)
        if (++debug_count, deeper(i->extra))
          v.push(i);
    }
    (void)debug_count;
  }
};

void run(auto const *const group, auto &&deeper) {
  std::remove_reference_t<decltype(*group)>::depth_first(group, deeper);
}

struct DeleteGroup {
  template <class G> void operator()(G *const g) const {
    G::depth_first_delete(g);
  }
};

inline DeleteGroup constexpr delete_group;

template <class E, class I>
auto tree(I const first, I const last, auto &&z) noexcept {
  struct {
    uint64_t mask = ~uint64_t{};
    void shift() { mask <<= 2; }
    void halt() { mask = {}; }
  } state;
  using G = Group<I, E>;
  using P = std::shared_ptr<G>;

  // Turn every particle into a group.
  auto lift = [first, last]() -> G * {
    if (first == last)
      // Degeneracy 0 (no particles).
      return {};
    // Let a and b be iterators to the particles exactly one particle apart
    // (single-particle range).
    I b = first, a = b++;
    auto *const g = new G{a, b};
    if (b == last)
      // Degeneracy 1 (one particle).
      return g;
    // Ordinary case (two+ particles).
    // Let g and h be pointers to groups.
    auto *root = new G{a, b}; // (minimize work).
    root->child = g;
    do {
      ++a, ++b;
      g = g->sibling = new G{a, b};
    } while (b != last);
    return root;
  };
  auto *g = lift();
  if (!g)
    return P{};
  if (!g->sibling) {
    auto *root = new G{first, last};
    root->child = g;
    return P{root, delete_group};
  }

  auto prefix = [&z, &state](auto &&a) { return z(a, state.mask); };
  auto complete = [&prefix, &state](auto *g) {
    assert(g && !g->sibling);
    auto *l = g->child;
    if (!l)
      return;
    auto a = prefix(*l->first);
    auto same = true;
    for (auto *m = l->sibling; m;) {
      if (auto b = prefix(*m->first); a != b) {
        if (l->sibling == m) {
          l->sibling = m->sibling;
          m->sibling = {};
          g = g->sibling = m;
          m = l->sibling;
        } else {
          assert(!g->sibling);
          g = g->sibling = new G{l, m};
          l = m;
          m = l->sibling;
        }
        a = b;
      } else {
        same = false;
        m = m->sibling;
      }
    }
    if (same)
      return state.halt();
    assert(!g->sibling);
    g->sibling = new G{l, {}};
    return;
  };

  auto *root = new G{*g};
  root->clear_relationships();
  root->child = g;
  while (state.mask) {
    state.shift(), complete(root);
    auto *j = root;
    root = new G{*root};
    root->clear_relationships();
    root->child = j;
  }

  return P{root, delete_group};
}

} // namespace detail

using detail::run;
using detail::tree;

} // namespace dyn::bh32

#endif // GRASS_BARNES_HUT_H
