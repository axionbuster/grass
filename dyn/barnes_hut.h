#ifndef GRASS_BARNES_HUT_H
#define GRASS_BARNES_HUT_H

#include <complex>
#include <cstdint>
#include <memory>
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

auto tree(auto first, auto last, auto &&z);
struct DeleteGroup;

template <class I, class E> class Group {
  friend auto tree(auto, auto, auto &&);
  friend struct DeleteGroup;

  I first, last;
  Group *child{}, *sibling{};

public:
  E extra;

  ~Group() = default;

private:
  Group(I first, I last) : first{first}, last{last}, extra{first, last} {}
  Group(Group const &g) : first{g.first}, last{g.last}, extra{g.extra} {}
  Group(Group &&g) noexcept
      : first{std::move(g.first)}, last{std::move(g.last)},
        extra{std::move(g.extra)} {}
  Group &operator=(Group const &g) {
    if (this != &g) {
      first = g.first, last = g.last, extra = g.extra;
    }
    return *this;
  }
  Group &operator=(Group &&g) noexcept {
    if (this != &g) {
      first = std::move(g.first), last = std::move(g.last),
      extra = std::move(g.extra);
    }
    return *this;
  }

  static void preorder(auto g, auto &&f) {
    if (!g)
      return;
    auto v = std::vector{g};
    while (!v.empty()) {
      auto h = v.back();
      v.pop_back();
      if (auto a = h->child) {
        v.push_back(a);
        for (auto b = a->sibling; b; b = b->sibling)
          v.push_back(b);
      }
      f(h);
    }
  }
};

struct DeleteGroup {
  void operator()(auto g) const {
    decltype(g)::preorder(g, [](auto h) { return delete h; });
  }
};

inline DeleteGroup constexpr deleteGroup;

template <class I, class E>
std::unique_ptr<Group<I, E>, DeleteGroup> tree(I const first, I const last,
                                               auto &&z) {
  struct {
    uint64_t mask = ~uint64_t{};
    void shift() { mask <<= 2; }
  } state;
  typedef Group<I, E> Group;

  // Turn every particle into a group.
  if (first == last)
    return {};
  I b = first, a = b++;
  auto *g = new Group{a, b};
  if (b == last)
    return {g, deleteGroup};
  auto *h = g;
  do {
    ++a, ++b;
    auto *i = new Group{a, b};
    h->sibling = i;
    h = i;
  } while (b != last);

  // Let's move out.
  h = new Group{*g};
  auto *i = h;
  h->child = g;
  state.shift();
  bool count{};
  auto prefix = [&z, &state](auto &&a) { return z(a, state.mask); };
  auto *m = g;
  auto c = prefix(*m);
  for (auto *n = m->sibling; n; n = n->sibling) {
    auto d = prefix(*n);
    if (c != d) {
      i->sibling = new Group{m, n};
      i = i->sibling;
      c = d;
      count = true;
    }
  }
  if (!count) {
    h->child = {};
    deleteGroup(h);
  }

  // FIXME: Put in loop; comment.
}

} // namespace detail

using detail::run;
using detail::tree;

} // namespace dyn::bh32

#endif // GRASS_BARNES_HUT_H
