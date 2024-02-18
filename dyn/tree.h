#ifndef GRASS_TREE_H
#define GRASS_TREE_H

#include <array>
#include <cassert>
#include <complex>
#include <memory>
#include <variant>

namespace dyn::tree32 {

class Tree;

/// @brief WIP.
class Node {
  friend class Tree;
  typedef std::array<std::unique_ptr<Node>, 4> Kids;
  typedef std::pair<size_t, size_t> Range;
  std::variant<Kids, Range> variant;
  std::complex<float> xy, ll;
  float m{}, s{};

  void place(size_t p, Kids &kids, int precision, auto get_xy, auto get_mass) {
    // Compute geometric center of this.
    auto center = ll + 0.5f * std::complex{s, s};
    auto pxy = get_xy(p);
    auto gx = pxy.real() > center.real();
    auto gy = pxy.imag() > center.imag();
    auto qu = int(gx) | int(gy) << 1;
    kids[qu]->rewrite(p, precision - 1, get_xy, get_mass);
  }

  void rewrite(size_t p, int precision, auto get_xy, auto get_mass) {
    if (std::holds_alternative<Range>(variant)) {
      auto &[b, e] = std::get<Range>(variant);
      if (b == e) {
        b = p, e = p + 1;
        xy = get_xy(p), m = get_mass(p);
      } else if (precision) {
        Kids kids;
        for (auto &&k : kids)
          // FIXME: Set center, side length
          k = std::make_unique<Node>();
        for (auto q = b; q != e; q++)
          place(q, kids, precision, get_xy, get_mass);
        variant = std::move(kids);
        // Try again with the same particle index (p).
        rewrite(p, precision, get_xy, get_mass);
      } else {
        // Extend list.
        e = p + 1;
        // FIXME: Assess potential accuracy problem.
        // Compute new mass (m) and center of mass (xy).
        auto pm = get_mass(p);
        xy *= m;
        xy += get_xy(p) * pm;
        m += pm;
        xy /= m;
      }
    } else {
      assert(precision);
      place(p, std::get<Kids>(variant), precision - 1, get_xy, get_mass);
    }
  }
};

class Tree {
  std::unique_ptr<Node> _root;
};

} // namespace dyn::tree32

#endif //GRASS_TREE_H
