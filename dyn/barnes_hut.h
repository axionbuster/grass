#ifndef GRASS_BARNES_HUT_H
#define GRASS_BARNES_HUT_H

#include <array>
#include <cassert>
#include <complex>
#include <cstdint>
#include <memory>
#include <optional>
#include <queue>
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
void run(auto &&, auto &&);
struct DeleteGroup;

template <class I, class E> class Group {
  template <class, class J> friend auto tree(J, J, auto &&) noexcept;
  friend void run(auto &&, auto &&);
  friend struct DeleteGroup;

  I first, last;
  Group *child{}, *sibling{};

public:
  E extra;

  ~Group() = default;

private:
  Group(I const first, I const last)
      : first{first}, last{last}, extra{first, last} {}
  Group(I const first, I const last, E const extra)
      : first{first}, last{last}, extra{extra} {}

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

void run(auto &&group, auto &&deeper) {
  std::remove_reference_t<decltype(*group)>::depth_first(&*group, deeper);
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
    [[nodiscard]] explicit operator bool() const { return mask; }
    void shift() { mask <<= 2; }
  } state;
  using G = Group<I, E>;
  using P = std::shared_ptr<G>;

  // Check for degeneracies (0 or 1 particle cases)
  if (first == last)
    return P{};
  else if (auto f = first; ++f == last)
    return P{new G{first, last}, delete_group};

  // Two or more particles.

  auto q = std::deque<G *>{};

  // First, turn every particle into a group.
  for (auto f = first; f != last; ++f) {
    auto g = f;
    ++g, q.push_back(new G{f, g});
  }
  state.shift();

  // Now, make layers.

  // Find the Z-code with the lowest few bits zeroed out.
  auto prefix = [&state, &z](auto &&p) { return z(p, state.mask); };

  while (state) {
    decltype(q) q2{};
    // Scan the queue (q) and then bring every subarray of groups with the same
    // Z-prefix under a common parent.
    assert(q.size());
    auto top = q.front();
    q.pop_front();
    class B {
      /// First particle.
      I first;

      /// Earliest group.
      G *group0;

      /// Latest group.
      G *group1;

    public:
      B(I f, G *g) : first{f}, group0{g}, group1{g} {}

      /// Admit a group.
      void admit(G *const g) { group1 = g; }

      /// Create a common parent group to all the included groups.
      [[nodiscard]] G *make() const {
        // Test: many groups or one group?
        if (group0 == group1) {
          // One group. Don't allocate; reuse.
          return group1;
        }
        // Many groups.
        auto h = new G{first, group1->last};
        h->child = group0;
        return h;
      }

      [[nodiscard]] I get_first() const { return first; }
    } build{top->first, top};
    auto zf = prefix(*build.get_first());
    for (auto &&g : q) {
      auto zg = prefix(*g->first);
      if (zf == zg) {
        // Same group. Merely update.
        build.admit(g);
      } else {
        // New group (new prefix).
        q2.push_back(build.make());
        // Update future new group.
        build = {g->first, g};
        // This new group will have this prefix.
        zf = zg;
      }
    }
    // Unconditional runoff: Handle it.
    q2.push_back(build.make());
    // Create or override siblings relationships in new layer (q2).
    assert(q2.size());
    for (typename decltype(q2)::size_type i = 0; i < q2.size() - 1; i++)
      q2[i]->sibling = q2[i + 1];
    // Recognize lower level as children.
    if (q2.front() != top)
      q2.front()->child = top;
    // Swap around.
    std::swap(q, q2);
    // Next level or stop.
    state.shift();
  }

  // Create and then set up root node. Return it.
  assert(q.size());
  auto root = new G{first, last};
  root->child = q.front();
  return P{root, delete_group};
}

} // namespace detail

using detail::run;
using detail::tree;

} // namespace dyn::bh32

#endif // GRASS_BARNES_HUT_H
