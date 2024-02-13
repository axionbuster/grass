#ifndef GRASS_TABLE_H
#define GRASS_TABLE_H

#include <vector>

#include <circle.h>
#include <kahan.h>
#include <newton.h>
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

  /// @brief Perform an integration step.
  /// @param dt Step size [units: T].
  void step(float dt) noexcept {
    auto copy{*this};
    for (size_t i = 0; i < size(); i++) {
      auto &&p = (*this)[i];
      // accel:
      // Supposing that particle p is located instead at the position xy below,
      // what is the acceleration experienced by p due to all the other
      // particles (q)?
      auto accel = [&](auto xy) {
        dyn::Kahan<std::complex<float>> a;
        // NOTE: the nested loop begins here.
        for (auto &&q : *this) {
          if (&p != &q) {
            dyn::Circle<float> cp{xy, p.radius}, cq = q.circle();
            a += gr.field(cp, cq, G * q.mass);
          }
        }
        return a();
      };
      // step(): Integrate.
      // (Calls accel() above a few times with slightly different xy.)
      auto step = Integrator{p.xy, p.v};
      step.step(dt, accel);
      auto &q = copy[i];
      q.xy = step.y0, q.v = step.y1;
    }
    *this = copy;
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
