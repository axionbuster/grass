#include <raylib.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <sstream>
#include <vector>

#include <circle.h>
#include <halton.h>

#include "Table.h"
#include "irhall.h"

using namespace phy;

static int do_main() {
  auto constexpr PARTICLES_LIMIT = 1'000;
  auto constexpr RADIUS = 0.05f;
  auto constexpr MASS = 1.0f;
  auto constexpr too_far = [](Particle &p) {
    return std::abs(p.xy) > 5'000.0f;
  };

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
    return Particle{{0}, {0}, mass, radius};
  };

  InitWindow(600, 600, "Basic 1,000 particle demo (click to add particles)");

  // The physical table (store particles, etc.); backup.
  Table table, table0;

  // In demo mode (default), show three particles in a figure-8 shape.
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
  table0 = table;

  // Flag for information to show.
  struct {
    // Show FPS.
    bool fps : 1 {};
    // Show the number of particles.
    bool n_particles : 1 {};
    // Show debug information about the camera.
    bool cam : 1 {};
    // Decide whether any flag is set.
    [[nodiscard]] constexpr bool any() const {
      return fps || n_particles || cam;
    }
    // Rotate to the next option.
    constexpr void next() {
      if (!fps)
        fps = n_particles = cam = true;
      else if (!cam)
        fps = n_particles = cam = false;
      else
        fps = n_particles = true, cam = false;
    }
  } show;

  // Interactive variables and flags.
  // (By itself, it returns whether the simulation is in interactive mode
  // [true] or demo mode [false]).
  struct {
    // If set, then the user intends to interact.
    bool on : 1 {};
    // If set, skip spawning in this frame when user asks to spawn a particle.
    bool spawned_last_frame : 1 {};
    // Do not advance the simulation in this frame.
    bool freeze : 1 {};
    // Target FPS.
    unsigned short target_fps{90};
    // Last time the simulation began or reset, seconds.
    double last_sec{GetTime()};
    constexpr explicit operator bool() const { return on; }
    [[nodiscard]] constexpr float target_dt() const {
      return 1.0f / float(target_fps);
    }
  } interactive;

  SetTargetFPS(interactive.target_fps);
  Camera2D cam{}, cam0{};
  cam.zoom = 100.0f; // [px] / [L] where L = world length unit; px = pixels.
  cam.offset = {float(GetScreenWidth()) / 2.0f,
                float(GetScreenHeight()) / 2.0f};
  cam0 = cam;

  while (!WindowShouldClose()) {
    // Apply variable time.
    // FIXME: The integrators Yoshida and Verlet currently do not support
    //  variable time; numerical errors may occur (but not programming errors).
    float dt = interactive ? GetFrameTime() : interactive.target_dt();

    // Reset the simulation (R).
    if (IsKeyPressed(KEY_R))
      goto reset_sim;

    // Reset if in demo for long enough.
    if (!interactive && GetTime() - interactive.last_sec >= 30.0f)
      goto reset_sim;

    // Particles too far from the origin will be removed.
    std::erase_if(table, too_far);

    // Toggle debug options (T).
    if (IsKeyPressed(KEY_T))
      show.next();

    // Freeze or unfreeze with space bar.
    if (IsKeyPressed(KEY_SPACE))
      interactive.freeze = !interactive.freeze;

    // Pan with the right mouse button.
    if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) {
      auto u = GetMouseDelta();
      auto v = std::complex<float>{u.x, u.y};
      auto w = std::complex<float>{cam.target.x, cam.target.y};
      auto x = w - v / cam.zoom;
      auto y = Vector2{x.real(), x.imag()};
      cam.target = y;
    }

    // Zoom with the wheel if given.
    // dwheel: usually -1.0f or 1.0f with computer mouse [but may be different].
    if (auto dwheel = GetMouseWheelMove()) {
      // dwheel != 0

      auto u = GetMousePosition();
      auto v = GetScreenToWorld2D(u, cam);

      cam.offset = u;
      cam.target = v;

      auto constexpr ZOOM_INCR = 10.0f;
      cam.zoom += dwheel * ZOOM_INCR;
      cam.zoom = std::max(cam.zoom, ZOOM_INCR);
      cam.zoom = std::min(cam.zoom, 20.0f * ZOOM_INCR);
    }

    // Spawn particles when asked.
    // Also, set interactive, and lower the gravitational constant.
    if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
      // Set interactive.
      interactive.on = true;

      // When interactive, lower the gravitational constant.
      table.G = 0.015625f;

      // Make sure the mouse is moving quickly (pixels per frame).
      // (Prevent cramping).
      auto constexpr FAST_ENOUGH = 4.0f;
      auto delta = GetMouseDelta();
      if (interactive.spawned_last_frame &&
          std::hypot(delta.x, delta.y) < FAST_ENOUGH) {
        // Resume a normal course of action.
        interactive.spawned_last_frame = false;
        goto simulate;
      }

      if (table.size() >= PARTICLES_LIMIT)
        table.erase(table.begin());

      auto m = GetMousePosition();
      auto n = GetScreenToWorld2D(m, cam);
      auto o = std::complex<float>(n.x, n.y);
      auto p = random_particle();
      // Set location
      p.xy = o;
      table.push_back(p);

      interactive.spawned_last_frame = true;
    } else {
      interactive.spawned_last_frame = false;
    }

  simulate:
    // Do the simulation!
    if (!interactive.freeze) {
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

    // Draw all particles in the frame.
    BeginMode2D(cam);
    for (auto &&p : table) {
      auto cp = p.circle();
      DrawCircleV({cp.real(), cp.imag()}, cp.radius, WHITE);
    }
    EndMode2D();

    // Compose text and show it.
    {
      // The standard library understands how to format a complex number, but,
      // understandably, knows nothing about Raylib's custom vector types.
      auto constexpr v2c = [](Vector2 v) {
        return std::complex<float>{v.x, v.y};
      };

      std::stringstream buf;
      if (!interactive)
        buf << "(Demo; click anywhere to add particles.)\n";
      if (!show.any())
        buf << "R to reset; T for debug";
      if (show.fps)
        buf << "FPS: " << GetFPS() << '\n';
      if (show.n_particles)
        buf << "N: " << table.size() << '\n';
      if (show.cam)
        buf << "Zoom: " << cam.zoom << "\nTarget: " << v2c(cam.target)
            << "\nOffset: " << v2c(cam.offset) << '\n';
      auto str = buf.str();
      DrawText(str.c_str(), 16, 16, 20, WHITE);
    }
    EndDrawing();
    continue;

  reset_sim:
    show = {};
    interactive = {};
    cam = cam0;
    table = table0;
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
