#ifndef GRASS_BARNES_HUT_H
#define GRASS_BARNES_HUT_H

#include <algorithm>
#include <array>
#include <complex>
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
  for (auto &&w : W)
    // Spread the bits out with zeros in between.
    // (ex) w = 0b1011 -> 0b01'00'01'01
    for (auto &&h : H)
      w = (w | (w << h.shift)) & h.mask;

  // Imaginary first.
  return W[0] | (W[1] << 1);
}
} // namespace detail

/// @brief Compute the Morton (Z) code of a complex number xy assuming a squared
/// grid by pre-multiplying the factor `Precision` to each component of xy.
/// If either component is not a finite number after the scaling, the result
/// is {} (no value).
template <uint32_t Precision = 512>
std::optional<uint64_t> morton(std::complex<float> xy) {
  xy *= float(Precision);
  if (std::abs(xy.real()) < float(INT32_MAX) &&
      std::abs(xy.imag()) < float(INT32_MAX)) {
    // unsigned.
    auto sgn = uint32_t(0x8000'0000ul);
    auto x = int32_t(xy.real()) ^ sgn; // promoted to unsigned.
    auto y = int32_t(xy.imag()) ^ sgn; // (ditto).
    return detail::interleave32(x, y);
  }
  return {};
}

/// @brief Group consecutive particles with the same Z-code.
/// @param begin Starting iterator for the particles.
/// @param end Past-the-end iterator for the particles, possibly equal to begin.
/// @param z Given what is like a const lvalue reference to a particle, compute
/// its z-code, perhaps masked (low bits flushed to zero).
/// @param grp Callback with two arguments, both iterators to the particles, to
/// be called when a range of particles is announced.
void group(auto begin, auto const end, auto &&z, auto &&grp) {
  while (begin != end) {
    auto [first, last] = std::ranges::equal_range(begin, end, z(*begin), {}, z);
    if (first != last)
      grp(first, begin = last);
    else
      ++begin;
  }
}

/// A consecutive group of particles.
/// @tparam E Extra data (type of field `data`). An E object may be constructed
/// by passing two iterators to particles ("begin" and "end) or by default
/// construction, and it can be updated using the `+=` operator with an E object
template <typename E, typename I> struct Group {
  I first, last;
  E data{};
  /// Create a group of particles using iterators to particles.
  Group(I first, I last) : first{first}, last{last}, data{first, last} {}
  Group(auto g_first, auto g_last) {
    first = g_first->first;
    last = g_last->last;
    while (g_first != g_last)
      data += g_first++->data;
  }
  [[nodiscard]] I begin() { return first; }
  [[nodiscard]] I end() { return last; }
  [[nodiscard]] I begin() const { return first; }
  [[nodiscard]] I end() const { return last; }
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

template <class E, class S> class View {
  /// Mask (begin with the finest detail, first).
  uint64_t mask = 0xffff'ffff'ffff'ffff;

  /// Iterate over particles.
  S s;

public:
  View(S s) : s(s) {}

  using Groups = std::vector<Group<E, decltype(s.begin())>>;

  /// Compute the groups at the current level of detail.
  /// @param z Take the particle and then compute the Morton code.
  /// @param prior The result of this function call for one finer level of
  /// detail.
  Groups groups(auto &&z, Groups const &prior = {}) const {
    if (!mask || s.begin() == s.end())
      return {};
    Groups novel{};
    auto m = mask;
    auto z_masked = [m, &z](auto &&p) -> std::optional<uint64_t> {
      if (auto w = z(p); w.has_value())
        return w.value() & m;
      else
        return {};
    };
    auto grp = [&novel](auto f, auto l) { novel.emplace_back(f, l); };
    if (prior.empty())
      return group(s.begin(), s.end(), z_masked, grp), novel;
    // One finer level of detail in `prior`.
    // Merge the groups, then, instead of recalculating everything.
    // Two-pointer solution: g and j are iterators to the groups in `prior`.
    auto g = prior.begin(), j = g;
    // Merger by having the same z-prefixes (a, b).
    auto a = z_masked(*g->first);
    while (++j != prior.end()) {
      auto b = z_masked(*g->first);
      // If same prefix, don't do anything specific.
      if (a != b) {
        // New prefix. Treat j as past-the-end group.
        novel.emplace_back(g, j);
        g = j;
        a = b;
      }
    }
    // Handle runoff.
    if (g != j)
      novel.emplace_back(g, j);
    return novel;
  }

  /// Make the level of detail one level coarser.
  void coarser() { mask <<= 2; }
};

template <class E, class S> auto levels(View<E, S> &view, auto &&z) {
  auto a = view.groups(z);
  auto l = std::vector<decltype(a)>{std::move(a)};
  do
    l.push_back(view.groups(z, l.back()));
  while (view.coarser(), !l.back().empty());
  l.pop_back();
  auto last = std::unique(l.begin(), l.end());
  l.erase(last, l.end());
  return l;
}

void run(auto &&levels, auto &&process) {
  if (levels.rbegin() != levels.rend()) {
    auto prior{*levels.rbegin()};
    auto begin = levels.rbegin();
    // Re-use allocated memory.
    decltype(prior) copy{};
    while (++begin != levels.rend()) {
      copy.clear();
      auto const &novel = *begin;
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
      // Process.
      std::erase_if(copy, process);
      // Break if no more.
      if (copy.empty())
        break;
      prior = copy;
    }
  }
}

} // namespace dyn::bh32

#endif // GRASS_BARNES_HUT_H
