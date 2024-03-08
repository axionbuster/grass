#include <circle.h>
#include <cmath>
#include <complex>
#include <halton.h>
#include <random>
#include <raylib.h>
#include <vector>

#include "Table.h"
#include "env.h"
#include "irhall.h"
#include "user.h"

using namespace phy;

/// Constants controlling the program.
struct Constants {
  /// Inclusive upper limit of the number of particles.
  size_t PARTICLES_LIMIT = 2'500;

  /// Uncorrected mean radius [L] and standard deviation [L].
  float LOG_MEAN_RADIUS = std::log(0.05f), LOG_STDEV_RADIUS = std::log(1.25f);

  /// Uncorrected mean mass [M] and standard deviation [M].
  float LOG_MEAN_MASS = std::log(1.0f), LOG_STDEV_MASS = std::log(1.0f);

  /// Squared distance for when a particle is too far [L^2].
  float SQ_DISTANCE_TOO_FAR = 5'000.0f * 5'000.0f;

  float G = 0.015625f;

  struct {
    bool galaxies : 1 {};
  } flags;

  /// Decide whether the position vector is too far.
  [[nodiscard]] constexpr bool too_far(std::complex<float> xy) const {
    return std::norm(xy) > SQ_DISTANCE_TOO_FAR;
  }

  [[nodiscard]] Particle random_particle(auto &&rng) const {
    using lognormal = std::lognormal_distribution<float>;
    lognormal m{LOG_MEAN_MASS, LOG_STDEV_MASS},
        r{LOG_MEAN_RADIUS, LOG_STDEV_RADIUS};
    // xy, v, m, r.
    return Particle{{}, {}, m(rng), r(rng)};
  }
};

template <typename... Args>
static constexpr Table<Args...> figure8(Constants constants) {
  Table table;

  // Make the mystical figure-8 shape below work at first.
  // (This G value is too large in most cases, so I will lower it once the
  // user starts interacting with the world, exiting the demo mode.)
  table.G = 1.0f;

  // Positions (c...), and velocities (v...)
  std::complex<float> c0{-0.97000436f, 0.24308753f},
      v0{0.4662036850f, 0.4323657300f}, v1{-0.93240737f, -0.86473146f};

  // Position, velocity, mass, radius
  table.emplace_back(c0, v0, constants.LOG_MEAN_MASS,
                     constants.LOG_MEAN_RADIUS);
  table.emplace_back(.0f, v1, constants.LOG_MEAN_MASS,
                     constants.LOG_MEAN_RADIUS);
  table.emplace_back(-c0, v0, constants.LOG_MEAN_MASS,
                     constants.LOG_MEAN_RADIUS);
  return table;
};

template <typename... Args>
static Table<Args...> galaxies(Constants constants) {
  using lognormal = std::lognormal_distribution<float>;
  using normal = std::normal_distribution<float>;
  using uniform = std::uniform_real_distribution<float>;

  auto const L = constants.PARTICLES_LIMIT;
  std::mt19937 rng{std::random_device{}()};
  struct {
    lognormal axes{-0.5f, 0.5f}, number{}, radial{};
    normal norm{};
    uniform angle{0.0f, 2.0f * std::numbers::pi_v<float>};
    [[nodiscard]] std::complex<float> normal_xy(std::mt19937 &rng) {
      return {norm(rng), norm(rng)};
    }
  } d{.number{lognormal{std::log(std::sqrt(float(L)))}}};

  Table table;
  table.G = constants.G;
  table.reserve(L);
  while (table.size() <= L) {
    auto const N = std::min(d.number(rng), float(L - table.size()));
    if (N <= 0.0f)
      break;
    auto first = table.size();
    for (size_t i = 0; i < size_t(N); i++)
      table.push_back(constants.random_particle(rng));
    auto ellipse = std::complex{d.axes(rng), d.axes(rng)};
    auto dot = [](auto a, auto b) -> std::complex<float> {
      return {a.real() * b.real(), a.imag() * b.imag()};
    };
    auto pan = d.normal_xy(rng) * 4.0f;
    auto spin = std::polar(4.0f, d.angle(rng));
    // Make an ellipse.
    for (auto i = first; i < table.size(); i++)
      table[i].xy = (dot(d.normal_xy(rng), ellipse) / 2.0f + pan) * spin;
  }
  return table;
}

static int do_main() {
  auto const constants = []() {
    Constants c;
    c.flags.galaxies = env::get("GRASS_GALAXIES").has_value();
    return c;
  }();

  std::mt19937 rng{std::random_device{}()};

  SetConfigFlags(FLAG_WINDOW_RESIZABLE);
  InitWindow(600, 600, "Basic 1,000 particle demo (click to add particles)");

  // Helper for the user interface.
  User user0;
  user0.control.demo = !constants.flags.galaxies;
  auto user{user0};

  auto const make_table = [&constants]() {
    return constants.flags.galaxies ? galaxies(constants) : figure8(constants);
  };

  // The physical table (store particles, etc.); backup.
  Table table = make_table();

  SetTargetFPS(user.control.target_fps);

  while (!WindowShouldClose()) {
    float dt = user.control.target_dt();

    // Reset the simulation (R) or if in demo for long enough.
    if (user.wants_reset() ||
        (user.control.demo && user.elapsed_sec() >= 30.0f))
      goto reset_sim;

    // Particles too far from the origin will be removed.
    std::erase_if(table,
                  [&constants](auto &&p) { return constants.too_far(p.xy); });

    // General interactions.
    user.rotate_debug_opts(), user.adjust_fly(), user.pan(), user.zoom();

    // Spawn particles when asked.
    // Also, clear user.control.demo, and lower the gravitational constant.
    if (auto xy = user.wants_spawn_particle(); xy.has_value()) {
      user.control.demo = false;

      // When user controls, reset the gravitational constant.
      table.G = constants.G;

      // Make sure the mouse is moving quickly (pixels per frame).
      // (Prevent cramping).
      auto constexpr FAST_ENOUGH = 4.0f;
      auto delta = GetMouseDelta();
      if (user.control.spawned_last_frame &&
          std::hypot(delta.x, delta.y) < FAST_ENOUGH) {
        // Resume a normal course of action.
        user.control.spawned_last_frame = false;
        goto simulate;
      }
      // Spawn a random particle at the mouse location.
      auto p = constants.random_particle(rng); // Mass and radius only.
      p.xy = xy.value();
      table.push_back(p);

      if (table.size() > constants.PARTICLES_LIMIT)
        table.erase(table.begin());

      user.control.spawned_last_frame = true;
    } else {
      user.control.spawned_last_frame = false;
    }

  simulate:
    // Do the simulation!
    if (user.control.fly) {
      table.step(dt);

      // Remove statistical bias in collision handling routine.
      // (See refresh_disk()'s comments for details.)
      table.refresh_disk();

      // Inspect for such things as NaN and Infinity.
      if (!table.good())
        // NaN or infinity somewhere. Reset the simulation.
        goto reset_sim;
    }

    BeginDrawing();
    ClearBackground(BLACK);

    // Draw all particles (p) visible in the window (w).
    BeginMode2D(user.cam);
    for (auto w = user.window(); auto &&p : table)
      if (auto c = p.circle(); dyn::intersect::disk_rectangle(c, w.ll, w.gg))
        user.particle(c);
    EndMode2D();

    // Compose text and show it.
    user.hud(table.size());
    EndDrawing();
    continue;

  reset_sim:
    user = user0;
    table = make_table();

    // (Do this or else hang.)
    BeginDrawing();
    EndDrawing();
  }
  CloseWindow();
  return 0;
}

#if defined(_WIN32)

#define U [[maybe_unused]]
int WinMain(U void **_0, U void **_1, U void **_2, U int _3) {
  return do_main();
}

#else

int main() { return do_main(); }

#endif
