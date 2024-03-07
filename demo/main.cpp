#include <circle.h>
#include <cmath>
#include <complex>
#include <halton.h>
#include <raylib.h>
#include <vector>

#include "Table.h"
#include "irhall.h"
#include "user.h"

using namespace phy;

static int do_main() {
  auto constexpr PARTICLES_LIMIT = 5'000;
  auto constexpr RADIUS = 0.05f;
  auto constexpr MASS = 1.0f;
  auto constexpr too_far = [](auto &&p) { return std::abs(p.xy) > 5'000.0f; };

  // Quasi-random number generator; uniform distribution on (0, 1).
  // - Quality not so important.
  // - Small state (only a single short int).
  // - Simple code (a single loop).
  // - Code is used in other places; possibility of inlining (as the compiler
  // judges is appropriate).
  dyn::Halton qrng;
  auto normal_variate = [&qrng]() {
    auto &q = qrng;
    // Approximate normal distribution using uniform random variates.
    return phy::irhnormal([&q]() { return q.x01(); });
  };
  // Make a particle at rest at the origin with a random mass and radius.
  auto random_particle = [&normal_variate]() {
    // Random mass and radius, but not other properties.
    auto mass = MASS * std::exp(normal_variate());
    auto radius = RADIUS * std::exp(2.0f * normal_variate());
    // xy, v, m, r.
    return Particle{{}, {}, mass, radius};
  };

  SetConfigFlags(FLAG_WINDOW_RESIZABLE);
  InitWindow(600, 600, "Basic 1,000 particle demo (click to add particles)");

  // Helper for the user interface.
  User user;

  // The physical table (store particles, etc.); backup.
  Table table;

  // In demo mode (default), user.show three particles in a figure-8 shape.
  {
    // Make the mystical figure-8 shape below work at first.
    // (This G value is too large in most cases, so I will lower it once the
    // user starts interacting with the world, exiting the demo mode.)
    table.G = 1.0f;

    // Positions (c...), and velocities (v...)
    std::complex<float> c0{-0.97000436f, 0.24308753f},
        v0{0.4662036850f, 0.4323657300f}, v1{-0.93240737f, -0.86473146f};

    // Position, velocity, mass, radius
    table.emplace_back(c0, v0, MASS, RADIUS);
    table.emplace_back(0.0f, v1, MASS, RADIUS);
    table.emplace_back(-c0, v0, MASS, RADIUS);
  }

  // Initial version of the table (backup).
  auto table0{table};

  SetTargetFPS(user.control.target_fps);

  while (!WindowShouldClose()) {
    float dt = user.control.target_dt();

    // Reset the simulation (R) or if in demo for long enough.
    if (user.wants_reset() ||
        (user.control.demo && user.elapsed_sec() >= 30.0f))
      goto reset_sim;

    // Particles too far from the origin will be removed.
    std::erase_if(table, too_far);

    // General interactions.
    user.rotate_debug_opts(), user.adjust_fly(), user.pan(), user.zoom();

    // Spawn particles when asked.
    // Also, clear user.control.demo, and lower the gravitational constant.
    if (auto xy = user.wants_spawn_particle(); xy.has_value()) {
      user.control.demo = false;

      // When user controls, lower the gravitational constant.
      table.G = 0.015625f;

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

      if (table.size() >= PARTICLES_LIMIT)
        table.erase(table.begin());

      // Spawn a random particle at the mouse location.
      auto p = random_particle(); // Mass and radius only.
      p.xy = xy.value();
      table.push_back(p);

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
    user = {};
    table = table0;

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
