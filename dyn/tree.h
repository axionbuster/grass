#ifndef GRASS_TREE_H
#define GRASS_TREE_H

#include <array>
#include <complex>
#include <memory>
#include <variant>

namespace dyn::tree {

/// @brief Expected interface for a particle.
template <typename P, typename F>
concept Particle = requires(P p, F f) {
  { p.get_xy() } -> std::convertible_to<std::complex<F>>;
  { p.get_mass() } -> std::convertible_to<F>;
};

template <std::floating_point F, Particle<F> P> class Tree;

class OutsideTree : public std::exception {
public:
  [[nodiscard]] char const *what() const override {
    return "particle seen outside tree boundaries";
  }
};

template <std::floating_point F, typename P> class Node {
  friend Tree<F, P>;

  /// @brief Variant active if this is an internal node.
  struct Internal {
    /// @brief Possibly null descendant nodes (pointers). "0b00" (imaginary
    /// bit is more significant) for the negative-negative quadrant, and
    /// "0b01" for the negative (imaginary)-positive (real) quadrant (etc.).
    /// So, the bits are in the "yx" order (if "y" is the imaginary).
    std::array<std::unique_ptr<Node>, 4> _descendants{};
  };

  /// @brief Relative location of the center of mass (this node is a square
  /// with two extreme corners at (-1, -1) at (1, 1), respectively.
  std::complex<F> _xy{};

  /// @brief Mass of the "impostor" particle.
  float _m{};

  /// @brief Hold either a (possibly null) pointer to the particle itself or the
  /// data structure in the case of an internal node. Do not manage the P object
  /// (if referencing). If the first variant, then either a terminal node (when
  /// the P pointer is not null) or an empty node (if null); If the second
  /// variant, an internal node.
  std::variant<P *, Internal> _x;

  enum struct Kind : signed char { empty, terminal, internal };

  /// Classify itself.
  [[nodiscard]] constexpr Kind kind() const {
    if (std::holds_alternative<P *>(_x)) {
      if (std::get<P *>(_x))
        return Kind::terminal;
      else
        return Kind::empty;
    }
    return Kind::internal;
  }

  Node() : _x(nullptr) {}
  Node(P *p, std::complex<F> xy, F mass) : _x(p), _xy(xy), _m(mass) {}

  void rewrite(P *p) {
    // p is valid.
    auto xy = p->get_xy();
    for (auto &&c : reinterpret_cast<F(&)[2]>(xy)) {
      if (c < F{-1} || c > F{1})
        throw OutsideTree{};
    }
    auto k = kind();
    if (k == Kind::empty) {
      _x = p;
      _xy = xy;
      _m = p->get_mass();
      return;
    }

    // Hard part begins...

    // Relative centers of the four quadrants.
    std::array<std::complex<F>, 4> constexpr centers{
        // Remember, imaginary component bit is more significant.
        // Bit 0 = negative; bit 1 = positive.
        // 0b00 ... yx
        {F{-0.5}, F{-0.5}}, // xy.
        // 0b01 ... yx
        {F{+0.5}, F{-0.5}}, // xy [tricky one].
        // 0b10 ... yx
        {F{-0.5}, F{+0.5}}, // xy.
        // 0b11 ... yx
        {F{+0.5}, F{0.5}}, // xy.
    };

    // Determine the quadrant of a point in the [-1,1]^2 square.
    auto constexpr quadrant = [](std::complex<F> xy) {
      return int(xy.real() < 0) | int(xy.imag() < 0) << 1;
    };

    if (k == Kind::terminal) {
      // Rewrite the proof (this) so that old particle (q) and new particle (p)
      // belong to the correct quadrant(s).
      auto *q = std::get<P *>(_x);
      // q is valid.
      //
    }
  }
};

/// @brief A quadtree centered about the origin.
template <std::floating_point F, Particle<F> P> class Tree {
  /// @brief (Possibly empty) reference to the root node.
  std::unique_ptr<Node<F, P>> _root;

  /// @brief Non-zero, positive, finite side length in world units [L].
  F _side;

  void rewrite(P &p) {
    auto xy = p.get_xy() / F{2} / _side + std::complex{F{0.5}, F{0.5}};
    if (!_root) {
      _root = std::make_unique<>(&p, xy, p.get_mass());
      return;
    }
  }
};

} // namespace dyn::tree

#endif // GRASS_TREE_H
