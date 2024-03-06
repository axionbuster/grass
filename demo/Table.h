#ifndef GRASS_TABLE_H
#define GRASS_TABLE_H

#include <algorithm>
#include <barnes_hut.h>
#include <cassert>
#include <circle.h>
#include <complex>
#include <cstdint>
#include <kahan.h>
#include <newton.h>
#include <optional>
#include <type_traits>
#include <vector>
#include <verlet.h>
#include <yoshida.h>

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
  dyn::Gravity<> gr;

  /// "Extra data" stored for a Barnes-Hut tree node. A circle.
  struct Physicals {
    /// Center [L].
    std::complex<float> xy;

    /// Radius [L] and mass [M].
    float radius{}, mass{};

    // Three required member functions (by `dyn::bh32::tree`):
    //  1. No-argument constructor.
    //  2. Constructor given a range of particles.
    //  3. Merger function using (+=).

    Physicals() = default;

    /// Given a range of particles (with an `xy` field), compute the quantities.
    template <class I> Physicals(I const first, I const last) {
      std::complex<double> xyd;
      for (auto i = first; i != last; ++i) {
        mass += i->mass;
        xyd += double(i->mass) * std::complex<double>{i->xy};
      }
      xy = std::complex<float>{xyd / double(mass)};
      for (auto i = first; i != last; ++i)
        radius = std::max(radius, std::abs(i->xy - xy));
    }

    /// Merge p's information.
    Physicals &operator+=(Physicals const &p) {
      if (this == &p)
        // Cannot violate const contract.
        return mass += mass, *this;
      // Compute the new average xy.
      auto sum = mass + p.mass;
      xy = mass / sum * xy + p.mass / sum * p.xy;
      // The rest.
      mass += p.mass;
      radius = std::max(radius, p.radius + std::abs(p.xy - xy));
      return *this;
    }
  };

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

    // Make a copy of this instance, and then set up a view.
    auto const copy{*this};
    // Apply bitwise AND with the mask (m) to the particle (p).
    auto morton_masked = [](auto &&p, auto m) -> std::optional<uint64_t> {
      if (auto z = p.morton; z.has_value())
        return z.value() & m;
      else
        return {};
    };
    // Compute the Barnes-Hut tree over the particles this has.
    auto const tree = bh::tree<Physicals>(begin(), end(), morton_masked);

    // Iterate over the particles, summing up their forces.
    for (size_t i = 0; i < size(); i++) {
      auto &&p = copy[i];
      // accel:
      // Supposing that particle p is located instead at the position xy below,
      // what is the acceleration experienced by p due to all the other
      // particles or approximations (g)?
      auto const tan_thr = tan_angle_threshold;
      auto const G = this->G;
      auto &&gr = this->gr;
      auto accel = [tan_thr, &gr, G, &p, &tree](auto xy) {
        std::complex<float> a;
        // These bool-coercing constants are due to
        // `dyn::bh32::Group<,>::depth_first`.
        enum { IGNORE = 0, DEEPER = 1 };
        // Look at the group (g) of particles.
        // Compute the acceleration and then exit the branch if far enough.
        // Otherwise, go deeper.
        auto deeper = [xy, tan_thr, &a, &gr, G, &p](auto &&g) {
          if (g.xy == xy)
            // Exactly equal? Must be self. Ignore.
            return IGNORE;
          auto distance = std::abs(g.xy - xy);
          if (distance < g.radius)
            return DEEPER;
          if (auto tan = g.radius / distance; tan_thr < tan)
            return DEEPER;
          // Treat group (g) as a point particle.
          // (Also, apply the gravitational constant [G] in a sneaky way
          // that won't stress the dynamic range in the middle of the
          // calculations.)
          a += gr.field(dyn::Circle{xy, p.radius}, dyn::Circle{g.xy, g.radius},
                        G * g.mass, distance);
          // (Enough detail.)
          return IGNORE;
        };
        // Compute the acceleration (a) on the tree.
        tree->depth_first(deeper);
        return a;
      };
      // step(): Integrate.
      // (Calls accel() above a few times with slightly different xy.)
      auto step = Integrator{p.xy, p.v};
      step.step(dt, accel);
      auto &q = (*this)[i];
      q.xy = step.y0, q.v = step.y1;
    }
  }

  /// @brief Refresh the "disk" used for parts of the calculation.
  void refresh_disk() noexcept { gr.refresh_disk(); }

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
