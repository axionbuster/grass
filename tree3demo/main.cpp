// Visualize the Morton code grouping algorithm.

#include <algorithm>
#include <array>
#include <barnes_hut.h>
#include <bit>
#include <cinttypes>
#include <complex>
#include <cstdio>
#include <morton.h>
#include <numbers>
#include <optional>
#include <random>
#include <raylib.h>
#include <vector>

#ifdef PI
// Raylib defines this macro.
#undef PI
#endif

struct Particle {
  std::complex<float> xy, v;
  static std::optional<uint64_t> morton(Particle const &n) {
    // 512 = "precision"
    return dyn::fixedmorton32<512>(n.xy);
  }
};

struct NodePhysical {
  std::complex<float> center;
  float radius{};
};

struct State {
  /// @brief Number of particles.
  static auto constexpr N = 1000;
  static auto constexpr PI = std::numbers::pi_v<float>;
  static auto constexpr ANGLES = std::array{
      PI / 24.0f, PI / 16.0f, PI / 12.0f, PI / 6.0f, PI / 3.0f,
  };
  std::vector<Particle> particles;
  std::vector<dyn::bh32::Node<NodePhysical, decltype(particles.begin())>> nodes;
  int64_t mask = 0xffff'ffff'ffff'0000;
  float angle_threshold = PI / 12.0f;
  bool fly{};
  static State fresh() {
    State s;
    // Make random particles.
    std::mt19937 rng(std::random_device{}());
    // mean, stddev
    std::normal_distribution<float> xy_dist{0.0f, 0.5f}, v_dist{0.0f, 0.25f};
    s.particles.reserve(N);
    for (int i = 0; i < N; i++)
      s.particles.emplace_back(std::complex{xy_dist(rng), xy_dist(rng)},
                               std::complex{v_dist(rng), v_dist(rng)});
    // [! highlight !] Make the circles.
    s.group();
    return s;
  }

  /// @brief Create the groups (the circles seen).
  void group() {
    nodes = {};
    // Sort by morton.
    std::ranges::sort(particles.begin(), particles.end(), {}, Particle::morton);
    // Compute groups/nodes.
    auto m = std::bit_cast<uint64_t>(mask);
    auto z_masked = [m](auto p) {
      auto xy = p.xy;
      auto z = dyn::fixedmorton32<512>(xy);
      if (z.has_value())
        return std::optional<uint64_t>{z.value() & m};
      else
        return std::optional<uint64_t>{};
    };
    auto with_node = [this](auto node) { nodes.push_back(node); };
    dyn::bh32::group<NodePhysical>(particles.begin(), particles.end(), z_masked,
                                   with_node);
    // Assign the extra data (center and radius).
    for (auto &&n : nodes) {
      std::complex<float> center;
      float radius{}, count{};
      for (auto &&p : n)
        // Welford's online mean algorithm; compute center.
        center += (p.xy - center) / ++count;
      for (auto &&p : n)
        // Compute radius.
        radius = std::max(radius, std::abs(p.xy - center));
      n.extra.center = center;
      n.extra.radius = radius;
    }
  }

  // Each quadrant requires two bits.

  void mask_left() {
    if (!(mask <<= 2))
      mask = 0xc000'0000'0000'0000;
  }

  void mask_right() { mask >>= 2; }
  void up_angle() { angle_next(ANGLES.begin(), ANGLES.end()); }
  void down_angle() { angle_next(ANGLES.rbegin(), ANGLES.rend()); }

private:
  State() = default;
  void angle_next(auto const begin, auto const end) {
    // Find the angle closest to `angle_threshold` in the list [begin, end).
    auto difference = [this](auto a) { return std::abs(a - angle_threshold); };
    auto closest = std::ranges::min_element(begin, end, {}, difference);
    if (closest != end && ++closest != end)
      angle_threshold = *closest;
  }
};

int do_main() {
  SetConfigFlags(FLAG_WINDOW_RESIZABLE);
  InitWindow(600, 600, "A");
  SetTargetFPS(60);
  auto s = State::fresh();
  Camera2D cam{};
  float shortest_zoom;
  {
    // Set camera.
    // (Do once at startup.)
    auto w = float(GetScreenWidth()), h = float(GetScreenHeight()),
         z = std::min(w, h);
    cam.offset = {w * 0.5f, h * 0.5f};
    cam.zoom = 0.25f * z;
    shortest_zoom = cam.zoom;
  }

  while (!WindowShouldClose()) {
    auto dt = 1.0f / std::max(float(GetFPS()), 40.0f);

    // Reset (R)
    if (IsKeyPressed(KEY_R)) {
      s = State::fresh();
      BeginDrawing();
      EndDrawing();
      continue;
    }

    // Shift mask (Left or right)
    if (IsKeyPressed(KEY_LEFT)) {
      s.mask_left();
      s.group();
    } else if (IsKeyPressed(KEY_RIGHT)) {
      s.mask_right();
      s.group();
    }

    // Change angle threshold (Up or down)
    if (IsKeyPressed(KEY_UP)) {
      s.up_angle();
    } else if (IsKeyPressed(KEY_DOWN)) {
      s.down_angle();
    }

    // Pan (right mouse)
    if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) {
      auto u = GetMouseDelta();
      auto v = std::complex{u.x, u.y};
      auto w = std::complex{cam.target.x, cam.target.y};
      auto x = w - v / cam.zoom;
      cam.target = {x.real(), x.imag()};
    }

    // Allow particles to fly or stop the particles (Space).
    if (IsKeyPressed(KEY_SPACE))
      s.fly = !s.fly;

    // Zoom (mouse wheel)
    if (auto wheel = GetMouseWheelMove()) {
      auto u = GetMousePosition();
      auto v = GetScreenToWorld2D(u, cam);
      cam.offset = u;
      cam.target = v;
      auto constexpr ZOOM_INCR = 5.0f;
      cam.zoom += wheel * ZOOM_INCR;
      cam.zoom = std::max(cam.zoom, 0.25f * shortest_zoom);
      cam.zoom = std::min(cam.zoom, 10.0f * shortest_zoom);
    }

    BeginDrawing();
    ClearBackground(BLACK);

    // Print top-left text.
    {
      auto mask = std::bit_cast<uint64_t>(s.mask);
      auto zeroes = std::countr_zero(mask);
      // 16: offset (y; px).
      // 24: 20 [= font size; px] + 4 [= padding; px].
      auto line_y = [](auto i) { return 16 + 24 * i; };
      auto i = 0;
      auto constexpr COLOR = WHITE;
      char text[128]{};
      std::snprintf(text, sizeof text, "Mask = 0x%" PRIx64 " (%d zeroes)", mask,
                    zeroes);
      DrawText(text, 16, line_y(i++), 20, COLOR);
      std::snprintf(text, sizeof text, "(%d nodes/groups)",
                    int(s.nodes.size()));
      DrawText(text, 16, line_y(i++), 20, COLOR);
      std::snprintf(text, sizeof text, "Angle threshold ~ %.1f deg",
                    s.angle_threshold * 180.0f / std::numbers::pi_v<float>);
      DrawText(text, 16, line_y(i++), 20, COLOR);
      DrawText("Left or right key to shift mask", 16, line_y(i++), 20, COLOR);
      DrawText("Up or down key to change angle threshold", 16, line_y(i++), 20,
               COLOR);
      DrawText("R to reset", 16, line_y(i++), 20, COLOR);
    }

    // Draw circles, particles, visualizations.
    BeginMode2D(cam);
    {
      // Mouse.
      auto [mx, my] = GetScreenToWorld2D(GetMousePosition(), cam);
      auto mouse = std::complex{mx, my};
      auto mouse_down = IsMouseButtonDown(MOUSE_BUTTON_LEFT);

      // 2 px.
      auto radius = 2.0f / cam.zoom;

      // Each node spans a range of particles.
      for (auto &&node : s.nodes) {
        // Draw the particles in the range.
        for (auto &&p : node) {
          auto &&xy = p.xy;
          DrawCircleV({xy.real(), xy.imag()}, radius, WHITE);
        }

        // Inspect the extra data (e).
        //  --> The circles.
        auto &&e = node.extra;
        if (e.radius) {
          auto center = Vector2{e.center.real(), e.center.imag()};
          auto distance = std::abs(mouse - e.center);
          if (mouse_down) {
            // Angle rejection visualization (left mouse).
            if (distance < e.radius) {
              DrawCircleLinesV(center, e.radius, GRAY);
            } else {
              // Mouse outside the circle.
              // Accept or reject? Approximate view angle by perpendicular
              // construction of the radius from the circle's center. Always
              // under-approximates true view angle (but is good).
              //
              // Nota bene:
              //
              // Let x, y, x0, y0, and m be complex numbers (m != 0). Let x
              // and y be variables. Then, the equation
              //
              //  y - y0 = m * (x - x0)
              //
              // models a general transformation between x and y involving
              // scaling (zoom), rotation (spin), and translation (pan).

              auto rc = e.center - mouse;
              auto rc_mag = std::abs(rc);
              auto rad_base = std::complex{0.0f, e.radius / rc_mag};
              auto perp_rad1 = rad_base * rc;
              auto u_perp_rad1 = perp_rad1 + e.center;
              auto perp_rad2 = -rad_base * rc;
              auto u_perp_rad2 = perp_rad2 + e.center;
              auto v_perp_rad1 =
                  Vector2{u_perp_rad1.real(), u_perp_rad1.imag()};
              auto v_perp_rad2 =
                  Vector2{u_perp_rad2.real(), u_perp_rad2.imag()};
              auto angle = 2.0f * std::atan2(e.radius, rc_mag);
              auto good = angle < s.angle_threshold;
              auto primary_color = good ? YELLOW : RED;
              DrawCircleLinesV(center, e.radius, primary_color);
              if (good) {
                // NaN -> !good -> this branch not hit.

                // Cast rays.
                auto line_color = Fade(primary_color, 0.5f);
                DrawLineV({mx, my}, v_perp_rad1, line_color);
                DrawLineV({mx, my}, v_perp_rad2, line_color);
              }
            }
          } else {
            // Ordinary view (no left mouse).
            DrawCircleLinesV(center, e.radius, WHITE);
            if (distance < e.radius) {
              // On mouse hover, fill circle with gray.
              auto fade = 0.75f * (1.0f - distance / e.radius);
              auto color = Fade(WHITE, fade);
              DrawCircleV(center, e.radius, color);
            }
          }
        }
      }
    }

    // If flight is enabled, let the particles evolve.
    if (s.fly) {
      for (auto &&p : s.particles) {
        auto &&xy = p.xy;
        xy += p.v * dt;
      }
      s.group();
    }
    EndMode2D();
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
