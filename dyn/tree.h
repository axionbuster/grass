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
  std::variant<Range, Kids> variant;
  std::complex<float> xy, ll;
  float m{}, s{};

  /// @brief Place the particle (p) in the quadrant set (kids) at the given
  /// level of precision. Call `rewrite` recursively with a lower precision.
  void place(size_t p, int precision, auto const &get_xy) {
    assert(precision);
    // Compute geometric center of this.
    auto center = ll + 0.5f * std::complex{s, s};
    auto pxy = get_xy(p);
    auto gx = pxy.real() > center.real();
    auto gy = pxy.imag() > center.imag();
    auto qu = int(gx) | int(gy) << 1;
    auto &kids = std::get<Kids>(variant);
    kids[qu]->rewrite(p, precision - 1, get_xy);
  }

  /// @brief Assuming that particles are in contiguous, consecutive ordering,
  /// insert the particle (p) at the given level of precision.
  void rewrite(size_t p, int const precision, auto const &get_xy) {
    // Action Plan.
    //
    // from  | stimulus       | to   | what to do
    // -----------------------------------
    // empty | precision == 0 | list | Place the first particle (p).
    // empty | precision != 0 | list | (Ditto.)
    // list  | precision == 0 | list | Extend the list with p.
    // list  | precision != 0 | kids | Break up, apply recursion.
    // kids  | precision == 0 | ---- | [!] Impossible or invalid.
    // kids  | precision != 0 | kids | Place p at the correct quadrant ("kid").

    if (std::holds_alternative<Kids>(variant))
      // "kids" case.
      return place(p, precision, get_xy);

    // "empty" or "list" case.
    auto &[b, e] = std::get<Range>(variant);
    if (b == e) {
      // "empty".
      // ==> Place the first particle (p).
      b = p, e = p + 1;
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
      // "ll": less-less corner of the square representing "this."
      kids[0b00]->ll = ll;
      kids[0b01]->ll = ll + C{half, 0.0f};
      kids[0b10]->ll = ll + C{0.0f, half};
      kids[0b11]->ll = ll + C{half, half};

      // Copy the index values of the references b and e, which become invalid.
      auto bb = b, ee = e;

      // Commit to new type.
      variant = std::move(kids);

      // Move all particles (q) in the list [b...e] to the quadrant set (kids).
      for (auto q = bb; q != ee; q++)
        place(q, precision, get_xy);

      // Try again with the same particle (p) and precision.
      rewrite(p, precision, get_xy);
    } else {
      // "list" + "precision == 0"
      // ==> Extend the list with p.
      // Append p. Assume that p == e (list must be contiguous).
      assert(p == e);
      e++;
      // (END).
    }
  }

  /// @brief Compute mass and center of mass of this node and its descendants.
  void weigh(auto const &get_xy, auto const &get_mass) {
    xy = {}, m = {};
    if (std::holds_alternative<Kids>(variant)) {
      auto &&kids = std::get<Kids>(variant);
      for (auto &&k : kids) {
        xy += k->xy * k->m;
        m += k->m;
      }
      // If no particles, let nothing happen.
      if (m != 0.0f)
        xy /= m;
      return;
    }
    auto &[b, e] = std::get<Range>(variant);
    for (auto p = b; p != e; p++) {
      auto pxy = get_xy(p);
      auto pm = get_mass(p);
      xy += pxy * pm;
      m += pm;
    }
    if (m != 0.0f)
      xy /= m;
  }

public:
  /// @brief Kind of a node (mutually exclusive).
  enum Kind { EMPTY, TERMINAL, INTERNAL };

  /// @brief Classify itself.
  [[nodiscard]] constexpr Kind kind() const {
    if (std::holds_alternative<Kids>(variant))
      return INTERNAL;
    auto &[b, e] = std::get<Range>(variant);
    return b == e ? EMPTY : TERMINAL;
  }

  /// @brief Return a node to the quadrant (be it "empty" --- but not
  /// unique_ptr empty; if internal).
  std::unique_ptr<Node> &quadrant(unsigned char q) {
    if (kind() != INTERNAL)
      throw;
    return std::get<Kids>(variant).at(q);
  }

  /// @brief Return the inclusive-exclusive range of particles (if not
  /// internal).
  [[nodiscard]] std::pair<size_t, size_t> particles() const {
    if (kind() == INTERNAL)
      throw;
    return std::get<Range>(variant);
  }
};

/// @brief A Barnes-Hut tree (WIP).
class Tree {
  /// @brief The root node (if it exists).
  std::unique_ptr<Node> _root;
  Tree(std::unique_ptr<Node> root) : _root(std::move(root)) {}

public:
  Tree() = default;
  Tree(Tree &&tree) noexcept : _root(std::move(tree._root)) {}
  Tree &operator=(Tree &&tree) noexcept {
    if (this == &tree)
      return *this;
    _root = std::move(tree._root);
    return *this;
  }

  /// @brief Construct a tree given a Morton-order (Z-order) sorted list of
  /// particles.
  static Tree from_z_ordered(int precision, size_t n, std::complex<float> ll,
                             float side, auto const &get_xy,
                             auto const &get_mass) {
    auto root = std::make_unique<Node>();
    root->ll = ll, root->s = side;
    for (size_t i = 0; i < n; i++) {
      root->rewrite(i, precision, get_xy);
    }
    root->weigh(get_xy, get_mass);
    return {std::move(root)};
  }
};

} // namespace dyn::tree32

#endif //GRASS_TREE_H
