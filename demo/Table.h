#ifndef GRASS_TABLE_H
#define GRASS_TABLE_H

#include <algorithm>
#include <barnes_hut.h>
#include <circle.h>
#include <complex>
#include <cstdint>
#include <newton.h>
#include <optional>
#include <type_traits>
#include <vector>
#include <verlet.h>

namespace phy {

struct Particle {
  /// @brief Kinematic properties.
  std::complex<float> xy, v;

  /// @brief Mass and radius.
  float mass = 1.0f, radius = 1.0f;

  /// @brief Latest Morton code (if any).
  std::optional<uint64_t> morton{};

  /// @brief Create a particle at rest at (0, 0) that has unit mass and radius.
  constexpr Particle() = default;

  /// @brief Construct a particle with the given position, velocity, mass, and
  /// radius, respectively.
  /// @param xy Position.
  /// @param v Velocity.
  /// @param m Mass, positive.
  /// @param r Radius, positive.
  constexpr Particle(std::complex<float> xy, std::complex<float> v, float m = 1,
                     float r = 1)
      : xy{xy}, v{v}, mass{m}, radius{r}, morton{} {}

  /// @brief Construct a circle that represents this particle.
  [[nodiscard]] constexpr dyn::Circle<float> circle() const {
    return {xy, radius};
  }
};

/// @brief A type of integrator accepted by Table.
template <typename I, typename F>
concept IntegratorType = requires(I i, std::complex<F> c) {
                           { I{c, c} } -> std::convertible_to<I>;
                           { i.y0 } -> std::convertible_to<std::complex<F>>;
                           { i.y1 } -> std::convertible_to<std::complex<F>>;
                           // Also, with a function f of type std::complex<F> ->
                           // std::complex<F>, and an F type value h, possible
                           // to advance internal state with the syntax:
                           //    i.step(h, f);
                           // [h ostensibly stands for "step size," F for
                           // "floating point," and f computes the second
                           // derivative from the zeroth derivative.]
                         };

/// @brief Store a vector of particles and integrate them using the provided
/// integrator type.
template <typename Integrator = dyn::Verlet<float>>
requires IntegratorType<Integrator, float>
class Table : public std::vector<Particle> {
  /// @brief Gravitational interaction.
  dyn::Gravity<> gravity;

  /// "Extra data" stored for a Barnes-Hut tree node. A circle.
  template <class I> struct Physicals {
    /// Center [L].
    std::complex<float> xy;

    /// Radius [L] and mass [M].
    float radius{}, mass{};

    /// First particle.
    I first;

    /// Has many particles?
    bool many{};

    // Three required member functions (by `dyn::bh32::Group`):
    //  1. No-argument constructor.
    //  2. Constructor given a range of particles.
    //  3. Merger function using (+=).

    Physicals() = default;

    /// Given a range of particles (with an `xy` field), compute the quantities.
    Physicals(I const first, I const last) : first{first} {
      std::complex<double> xyd;
      unsigned char count{}; // (Only used to decide many particles vs. single.)
      for (auto i = first; i != last; ++i) {
        mass += i->mass;
        xyd += double(i->mass) * std::complex<double>{i->xy};
        using C = unsigned char;
        count = std::min(C(count + C(1)), C(2));
      }
      assert(count);
      many = count > 1;
      xy = std::complex<float>{xyd / double(mass)};
      for (auto i = first; i != last; ++i)
        radius = std::max(radius, i->radius + std::abs(i->xy - xy));
    }

    /// Merge p's information.
    Physicals &operator+=(Physicals const &p) {
      if (this == &p)
        // Cannot violate const contract.
        return mass *= 2.0f, *this;
      // Compute the new average xy.
      auto sum = mass + p.mass;
      xy = mass / sum * xy + p.mass / sum * p.xy;
      // The rest.
      mass += p.mass;
      radius = std::max(radius, p.radius + std::abs(p.xy - xy));
      // No need to update `first`:
      // Assume that mergers come "in order."
      return *this;
    }

    /// Non-required method to create a `Circle` instance.
    [[nodiscard]] dyn::Circle<float> circle() const { return {xy, radius}; }
  };

  /// @brief Given a Barnes-Hut tree and a circle that represents a particle,
  /// compute the acceleration onto the particle due to the data in the tree.
  std::complex<float> accelerate(auto &&tree, dyn::Circle<> circle, auto i) {
    std::complex<float> a{};
    tree->depth_first([this, circle, i, &a](auto &&group) {
      auto const TRUNCATE = false;
      auto const square = [](auto x) { return x * x; };
      if (!group.many && group.first == i)
        // Exclude self-interactions.
        return TRUNCATE;
      auto norm = std::norm(group.xy - circle), rsq = square(group.radius);
      // If a non-singular group either:
      //  - contains the center of `circle` inside said group's circle, or
      //  - the (underapproximated) view angle is too wide, then
      // resolve more detail.
      if (group.many &&
          (norm < rsq || square(tan_angle_threshold) < rsq / norm))
        return !TRUNCATE;
      // Compute the acceleration due to the group.
      // Also, insert the value of G, the universal gravitational constant, in a
      // way that doesn't stress the single-precision dynamic range.
      a += gravity.field(circle, group.circle(), G * group.mass,
                         std::sqrt(norm));
      return TRUNCATE;
    });
    return a;
  }

public:
  /// @brief Universal gravitational constant [LLL/M/T/T]. Modify freely.
  float G{1.0f};
  float tan_angle_threshold{0.12278456f}; // tan(7 deg)

  /// @brief Perform an integration step.
  /// @param dt Step size [units: T].
  void step(float dt) noexcept {
    namespace bh = dyn::bh32;

    // Compute the Morton code of all particles.
    for (auto &&p : *this)
      p.morton = bh::morton(p.xy);
    // Sort the particles in Z-order.
    std::ranges::stable_sort(begin(), end(), {},
                             [](auto &&p) { return p.morton; });

    // Apply bitwise AND with the mask (m) to the particle (p).
    auto morton_masked = [](auto &&p, auto m) -> std::optional<uint64_t> {
      if (auto z = p.morton; z.has_value())
        return z.value() & m;
      else
        return {};
    };
    // Compute the Barnes-Hut tree over the particles this has.
    auto const tree =
        bh::tree<Physicals<decltype(begin())>>(begin(), end(), morton_masked);

    // Iterate over the particles, summing up their forces.
    for (auto i = begin(); i != end(); ++i) {
      auto &&p = *i;
      // Supposing that particle p is located instead at the position xy below,
      // what is the acceleration experienced by p due to all the other
      // particles or approximations (g)?
      auto ig = Integrator{p.xy, p.v};
      ig.step(dt, [this, &tree, &p, &i](auto xy) {
        return this->accelerate(tree, {xy, p.radius}, i);
      });
      p.xy = ig.y0, p.v = ig.y1;
    }
  }

  /// @brief Refresh the "disk" used for parts of the calculation.
  void refresh_disk() noexcept { gravity.refresh_disk(); }

  /// @brief Test whether the simulation is in "good state."
  bool good() noexcept {
    auto constexpr finite = [](std::complex<float> f) {
      return std::isfinite(f.real()) && std::isfinite(f.imag());
    };
    for (auto &&p : *this)
      if (!finite(p.xy) || !finite(p.v))
        return false;
    return true;
  }
};

} // namespace phy

#endif // GRASS_TABLE_H
