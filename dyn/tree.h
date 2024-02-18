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

  /// @brief Place the particle (p) in the quadrant set (kids) at the given
  /// level of precision. Call `rewrite` recursively with a lower precision.
  void place(size_t p, Kids &kids, int precision, auto get_xy, auto get_mass) {
    assert(precision > 0);
    // Compute geometric center of this.
    auto center = ll + 0.5f * std::complex{s, s};
    auto pxy = get_xy(p);
    auto gx = pxy.real() > center.real();
    auto gy = pxy.imag() > center.imag();
    auto qu = int(gx) | int(gy) << 1;
    kids[qu]->rewrite(p, precision - 1, get_xy, get_mass);
  }

  /// @brief Assuming that particles are in contiguous, consecutive ordering,
  /// insert the particle (p) at the given level of precision (non-negative).
  void rewrite(size_t p, int const precision, auto get_xy, auto get_mass) {
    // Action Plan.
    //
    // type  | stimulus       | what to do
    // -----------------------------------
    // empty | precision == 0 | Place the first particle (p).
    // empty | precision != 0 | (Ditto.)
    // list  | precision == 0 | Extend the list with p.
    // list  | precision != 0 | Break up into quadrants, apply recursion.
    // quad  | precision == 0 | [!] Impossible or invalid.
    // quad  | precision != 0 | Place p at the correct quadrant.

    if (std::holds_alternative<Kids>(variant))
      // "quad" case.
      return place(p, std::get<Kids>(variant), precision, get_xy, get_mass);

    // "empty" or "list" case.
    auto &[b, e] = std::get<Range>(variant);
    if (b == e) {
      // "empty".
      // ==> Place the first particle (p).
      b = p, e = p + 1;
      xy = get_xy(p), m = get_mass(p);
    } else if (precision) {
      // "list" + "precision != 0"
      // ==> Break up into quadrants, apply recursion.

      // Allocate and initialize the quadrants ("kids").
      // - Set up the side length (s) and less-less corner (ll).
      Kids kids;
      auto half = s * 0.5f;
      for (auto &&k : kids)
        k = std::make_unique<Node>(), k->s = half;
      typedef std::complex<float> C;
      // significant bit: y (1 if greater than center y).
      // less significant bit: x (1 if greater than center x).
      kids[0b00]->ll = ll;
      kids[0b01]->ll = ll + C{half, 0.0f};
      kids[0b10]->ll = ll + C{0.0f, half};
      kids[0b11]->ll = ll + C{half, half};

      // Move all particles (q) in the list [b...e] to the quadrant set (kids).
      for (auto q = b; q != e; q++)
        place(q, kids, precision, get_xy, get_mass);

      // Commit to new type.
      variant = std::move(kids);

      // Try again with the same particle (p) and precision.
      rewrite(p, precision, get_xy, get_mass);
    } else {
      // "list" + "precision == 0"
      // ==> Extend the list with p.

      // Append p. Assume that p == e + 1.
      assert(p == e + 1);
      e++;
      // FIXME: Assess potential accuracy problem.
      // Compute new mass (m) and center of mass (xy).
      auto pm = get_mass(p);
      xy *= m;
      xy += get_xy(p) * pm;
      m += pm;
      xy /= m;
    }
  }
};

class Tree {
  std::unique_ptr<Node> _root;
};

} // namespace dyn::tree32

#endif //GRASS_TREE_H
