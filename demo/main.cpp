#include <raylib.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <sstream>
#include <vector>

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

  // Quasi-random number generator.
  // - Quality not so important.
  // - Possibly reduce binary size by reusing code.
  // - Possibly help with instruction cache by having simple code.
  // - Definitely small state (1 short int vs. compared to, say, mt19937).
  dyn::Halton qrng;
  auto normal_variate = [&qrng]() {
    auto &q = qrng;
    // Approximate normal distribution using uniform random variates.
    return phy::irhnormal([&q]() { return q.x01(); });
  };
  // Make a particle at rest at the origin with a random mass and radius.
  auto random_particle = [&normal_variate]() {
    // Random mass and radius, but not other properties.
    auto mass = std::exp(normal_variate());
    auto radius = std::exp(2.0f * normal_variate() - 3.0f);
    // xy, v, m, r.
    return Particle{{0}, {0}, mass, radius};
  };

  InitWindow(600, 600, "Basic 1,000 particle demo (click to add particles)");

  // The physical table (store particles, etc.); backup.
  Table table, table0;

  // Make the mystical figure-8 shape below work at first.
  // (This G value is too large in most cases, so lower it once the user starts
  // interacting with the world.)
  table.G = 1.0f;

  // Get started with three particles in a figure-8 shape.
  {
    // Position, velocity.
    std::complex<float> _c0{-0.97000436f, 0.24308753f},
        _v0{0.4662036850f, 0.4323657300f}, _v1{-0.93240737f, -0.86473146f};

    // Position, velocity, mass, radius
    table.emplace_back(_c0, _v0, MASS, RADIUS);
    table.emplace_back(0.0f, _v1, MASS, RADIUS);
    table.emplace_back(-_c0, _v0, MASS, RADIUS);
  }

  // Initial version of the table (backup).
  table0 = table;

  // Information to show or not.
  struct {
    bool fps : 1 {};
    bool n_particles : 1 {};
    bool cam : 1 {};
    [[nodiscard]] constexpr bool any() const {
      return fps || n_particles || cam;
    }
    constexpr void next() {
      if (!fps)
        fps = n_particles = cam = true;
      else if (!cam)
        fps = n_particles = cam = false;
      else
        fps = n_particles = true, cam = false;
    }
  } show;

  // If set, then the user intends to interact.
  struct {
    bool on : 1 {};
    bool spawned_last_frame : 1 {};
    bool other_user_manip : 1 {};
    bool freeze : 1 {};
    unsigned short target_fps{90};
    double last_sec{GetTime()};
    constexpr explicit operator bool() const { return on; }
    [[nodiscard]] constexpr float target_dt() const {
      return 1.0f / float(target_fps);
    }
  } interactive;

  SetTargetFPS(interactive.target_fps);
  Camera2D cam{}, cam0{};
  cam.zoom = 100.0f;
  cam.offset = {float(GetScreenWidth()) / 2.0f,
                float(GetScreenHeight()) / 2.0f};
  cam0 = cam;

  while (!WindowShouldClose()) {
    float dt = interactive ? GetFrameTime() : interactive.target_dt();

    // Reset the simulation (R) if asked.
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

    // Zoom with the wheel.
    if (float wheel = GetMouseWheelMove()) {
      auto u = GetMousePosition();
      auto v = GetScreenToWorld2D(u, cam);

      cam.offset = u;
      cam.target = v;

      auto constexpr ZOOM_INCR = 10.0f;
      cam.zoom += wheel * ZOOM_INCR;
      cam.zoom = std::max(cam.zoom, ZOOM_INCR);
      cam.zoom = std::min(cam.zoom, 20 * ZOOM_INCR);
    }

    // Spawn particles when asked.
    // Also, set interactive, and lower the gravitational constant.
    if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
      // Set interactive.
      interactive.on = true;

      // The user is manipulating it right now.
      interactive.other_user_manip = true;

      // When interactive, lower the gravitational constant.
      table.G = 0.015625f;

      // Make sure the mouse is moving quickly (pixels per frame).
      // (Prevent cramping).
      auto constexpr FAST_ENOUGH = 4.0;
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

      // The user has stopped manipulating.
      interactive.other_user_manip = false;
    }

  simulate:
    // Do the simulation!
    if (!interactive.freeze)
      table.step(dt);

    // table.refresh_disk();

    // Don't be moving the particles while it's being manipulated.
    // if (!interactive.spawned_last_frame && !interactive.other_user_manip)
      // table.center();

    // Inspect floating point exceptions.
    if (!table.good())
      // NaN or infinity somewhere. Reset the simulation.
      goto reset_sim;

    BeginDrawing();
    ClearBackground(BLACK);

    // Draw all particles in the frame.
    BeginMode2D(cam);
    for (auto &&p : table) {
      auto cp = p.circle();
      Vector2 xy{cp.real(), cp.imag()};
      // Should I issue the draw call or not? Tell it here:
      // co stores screen dimensions calculated earlier.
      auto co = cam.offset;
      // x, y, width, then height.
      Rectangle screen{-co.x / 2, -co.y / 2, co.x, co.y};
      if (CheckCollisionCircleRec(xy, cp.radius, screen))
        DrawCircleV(xy, cp.radius, WHITE);
    }
    EndMode2D();

    // Compose text and show it.
    {
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

int WinMain(void **_0, void **_1, void **_2, int _3) { return do_main(); }

#else

int main() { return do_main(); }

#endif
