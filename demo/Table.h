#ifndef GRASS_TABLE_H
#define GRASS_TABLE_H

#include <algorithm>
#include <barnes_hut.h>
#include <cassert>
#include <circle.h>
#include <kahan.h>
#include <newton.h>
#include <optional>
#include <vector>
#include <verlet.h>
#include <yoshida.h>

namespace phy {

struct Particle {
  /// @brief Kinematic properties.
  std::complex<float> xy, v;

  /// @brief Mass and radius.
  float mass = 1.0f, radius = 1.0f;

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
      : xy{xy}, v{v}, mass{m}, radius{r} {}

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
  // Also, with a function f of type std::complex<F> -> std::complex<F>,
  // and an F type value h,
  // possible to advance internal state with the syntax:
  //    i.step(h, f);
  // [h ostensibly stands for "step size," F for "floating point,"
  // and f computes the second derivative from the zeroth derivative.]
};

namespace detail {

struct W {
  std::complex<float> xy;
  float mass{}, radius{};
  W() = default;
  W(auto first, auto last) {
    for (auto p = first; p != last; ++p)
      mass += p->mass;
    if (auto p = first; p != last && ++p == last) {
      xy = first->xy;
      mass = first->mass;
      radius = first->radius;
      return;
    }
    std::complex<double> xyd;
    for (auto p = first; p != last; ++p) {
      radius = std::max(radius, p->radius + std::abs(p->xy - xy));
      xyd += double(p->mass) * std::complex<double>{p->xy};
    }
    xy = std::complex<float>{xyd};
  }
  W &operator+=(W const &w) {
    if (this == &w) {
      mass *= 2.0f;
      return *this;
    }
    auto xyd = double(mass) * std::complex<double>{xy} +
               double(w.mass) * std::complex<double>{w.xy};
    auto m = double(mass) + double(w.mass);
    xy = std::complex<float>{xyd / m};
    mass = float(m);
    radius = std::max(radius, w.radius + std::abs(w.xy - xy));
    return *this;
  }
};

} // namespace detail

/// @brief Store a vector of particles and integrate them using the provided
/// integrator type.
template <typename Integrator = dyn::Yoshida<float>>
requires IntegratorType<Integrator, float>
class Table : public std::vector<Particle> {
  template <typename F = double> using C = std::complex<F>;

  /// @brief Gravitational interaction.
  dyn::Gravity<> gr;

public:
  /// @brief Universal gravitational constant [LLL/M/T/T]. Modify freely.
  float G{1.0f};
  float tan_angle_threshold{0.0874887f}; // tan(5 deg)

  /// @brief Perform an integration step.
  /// @param dt Step size [units: T].
  void step(float dt) noexcept {
    namespace bh = dyn::bh32;

    // Sort the particles in Z-order.
    auto morton = [](auto &&p) { return bh::morton(p.xy); };
    auto morton_masked = [&morton](auto &&p,
                                   auto m) -> std::optional<uint64_t> {
      if (auto z = morton(p); z.has_value())
        return z.value() & m;
      else
        return {};
    };
    std::ranges::sort(begin(), end(), {}, morton);

    // Make a copy of this instance, and then set up a view.
    auto const copy{*this};
    bh::View<Table const &> view{copy};
    auto const levels = bh::levels<detail::W>(view, morton_masked);

    // Iterate over the particles, summing up their forces.
    for (size_t i = 0; i < size(); i++) {
      auto &&p = copy[i];
      // accel:
      // Supposing that particle p is located instead at the position xy below,
      // what is the acceleration experienced by p due to all the other
      // particles or approximations (g)?
      auto accel = [&](auto xy) {
        dyn::Kahan<std::complex<float>> a;
        enum { KEEP = 0, ERASE = 1 };
        auto process = [&](auto &&g) {
          auto gxy = g.data.xy;
          auto gra = g.data.radius;
          auto gm = g.data.mass;
          auto dist = std::abs(gxy - xy);
          if (g.single() && gxy == xy)
            return ERASE;
          if (dist < gra)
            return KEEP;
          if (auto tan = gra / dist; tan_angle_threshold < tan)
            return KEEP;
          auto cp = dyn::Circle{xy, p.radius}, cg = dyn::Circle{gxy, gra};
          a += gr.field(cp, cg, G * gm);
          return ERASE;
        };
        bh::run(levels, process);
        return a();
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
