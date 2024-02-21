// Visualize the Morton code grouping algorithm.

#include <algorithm>
#include <array>
#include <barnes_hut.h>
#include <bit>
#include <cinttypes>
#include <cmath>
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

/// @brief A node in the Barnes-Hut tree at a given depth (span zero or more
/// particles).
/// @tparam E Any extra data (default-constructible).
/// @tparam I Iterator type.
template <std::default_initializable E, typename I> struct Node {
  I first, last;
  E extra;
  I begin() { return first; }
  I begin() const { return first; }
  I end() { return last; }
  I end() const { return last; }
};

struct State {
  /// @brief Number of particles.
  static auto constexpr N = 1000;
  static auto constexpr PI = std::numbers::pi_v<float>;
  static auto constexpr ANGLES = std::array{
      PI / 24.0f, PI / 16.0f, PI / 12.0f, PI / 6.0f, PI / 3.0f,
  };
  std::vector<Particle> particles;
  std::vector<Node<NodePhysical, decltype(particles.begin())>> nodes;
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
    auto with_node = [this](auto first, auto last) {
      nodes.push_back({first, last, {}});
    };
    dyn::bh32::group(particles.begin(), particles.end(), z_masked, with_node);
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

    // Draw circles, particles, visualizations.

    // Count "good" nodes.
    struct {
      int good_nodes{};
      bool casting{};
    } ray_report;

    BeginMode2D(cam);
    {
      // Mouse.
      auto [mx, my] = GetScreenToWorld2D(GetMousePosition(), cam);
      auto mouse = std::complex{mx, my};
      auto mouse_down = IsMouseButtonDown(MOUSE_BUTTON_LEFT);
      ray_report.casting = mouse_down;

      // The tangent of the angle threshold.
      auto tan_angle_threshold = std::tan(s.angle_threshold);

      // Each node spans a range of particles.
      for (auto &&node : s.nodes) {
        // Draw the particles in the range.
        for (auto &&p : node) {
          auto &&xy = p.xy;
          auto radius = 2.0f / cam.zoom; // 2 px.
          DrawCircleV({xy.real(), xy.imag()}, radius, WHITE);
        }

        // Inspect the extra data (e).
        //  --> The circles.
        auto &&e = node.extra;
        if (auto radius = e.radius) {
          auto center = e.center;
          auto v_center = Vector2{center.real(), center.imag()};
          auto displacement = center - mouse;
          auto distance = std::abs(displacement);
          if (mouse_down) {
            // Angle rejection visualization (left mouse).
            if (distance < radius) {
              DrawCircleLinesV(v_center, radius, GRAY);
            } else {
              // Mouse outside the circle.
              //
              // Accept or reject? Approximate view angle by construction of
              // the endpoints of two radii perpendicular to the line of sight
              // from the mouse to the circle's center. Always
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

              // See the definition of the tangent:
              //  Suppose a line and a unit circle intersect exactly once. The
              //  line would therefore be called a tangent line of the circle.
              //  Now, exactly one radius of the circle intersects the tangent
              //  line. Then, draw any diameter of the circle and extend it on
              //  both sides to infinity. Suppose the extended line of the
              //  diameter intersects the tangent line. Then, a right triangle
              //  is formed between the radius, the tangent line, and the
              //  diameter line. The length of the triangle's side that lies on
              //  the tangent line is called the tangent (tan) of the angle (A)
              //  between the radius and the diameter line. (Let's call A the
              //  "central angle" from now on because A is formed at the center
              //  of the circle.) Notice: `tan A` is called "tangent" because it
              //  is tangent to the circle! Similarly, the diameter, when
              //  extended on both sides, would intersect the circle twice: it
              //  is a secant line. Among the sides of the right triangle,
              //  excluding the radius, only the hypotenuse is secant to the
              //  circle (when extended; this fact is best understood by drawing
              //  the figure). Hence, the hypotenuse of this triangle is called
              //  the secant (sec) of the central angle (A). Noticing that this
              //  right triangle has "1" and "tan A" as the legs and "sec A" as
              //  the hypoteneuse, by an elementary application of the
              //  Pythagorean theorem, it is seen that:
              //
              //    1 + tan^2 A = sec^2 A.
              //
              //  Hope this helps.
              auto tangent = e.radius / distance;
              auto c_tangent = std::complex{1.0f, tangent};

              auto rotate = displacement;
              auto pan = mouse;
              // Construct the radial endpoints perpendicular to displacement.
              auto r0 = c_tangent * rotate + pan;
              auto r1 = std::conj(c_tangent) * rotate + pan;
              auto good = tangent < tan_angle_threshold;
              auto primary_color = good ? YELLOW : RED;
              DrawCircleLinesV(v_center, radius, primary_color);
              if (good) {
                // NaN -> !good -> this branch not hit.
                ray_report.good_nodes++;
                // Cast rays.
                auto line_color = Fade(primary_color, 0.5f);
                DrawLineV({mx, my}, {r0.real(), r0.imag()}, line_color);
                DrawLineV({mx, my}, {r1.real(), r1.imag()}, line_color);
              }
            }
          } else {
            // Ordinary view (no left mouse).
            DrawCircleLinesV(v_center, e.radius, WHITE);
            if (distance < e.radius) {
              // On mouse hover, fill circle with gray.
              auto fade = 0.75f * (1.0f - distance / e.radius);
              auto color = Fade(WHITE, fade);
              DrawCircleV(v_center, e.radius, color);
            }
          }
        }
      }
    }
    EndMode2D();

    // Print top-left text.
    {
      auto mask = std::bit_cast<uint64_t>(s.mask);
      auto zeroes = std::countr_zero(mask);
      // 16: offset (y; px).
      // 24: 20 [= font size; px] + 4 [= padding; px].
      auto y = [](auto i) { return 16 + 24 * i; };
      // Line number (start at 0).
      auto i = 0;
      auto constexpr COLOR = WHITE;
      char text[128]{};
      // Count plural nodes (nodes with two or more particles).
      auto plural = 0;
      for (auto &&n : s.nodes)
        plural += int(n.last != n.first && n.last != n.first + 1);
      std::snprintf(text, sizeof text, "Mask = 0x%" PRIx64 " (%d zeroes)", mask,
                    zeroes);
      DrawText(text, 16, y(i++), 20, COLOR);
      if (ray_report.casting)
        std::snprintf(text, sizeof text, "(accept : reject = %d : %d)",
                      ray_report.good_nodes, plural - ray_report.good_nodes);
      else
        std::snprintf(text, sizeof text, "(%d nodes [plural %d nodes])",
                      int(s.nodes.size()), plural);
      DrawText(text, 16, y(i++), 20, COLOR);
      std::snprintf(text, sizeof text, "Angle threshold ~ %.1f deg",
                    s.angle_threshold * 180.0f / std::numbers::pi_v<float>);
      DrawText(text, 16, y(i++), 20, COLOR);
      DrawText("Left or right key to shift mask", 16, y(i++), 20, COLOR);
      DrawText("Up or down key to change angle threshold", 16, y(i++), 20,
               COLOR);
      DrawText("R to reset", 16, y(i++), 20, COLOR);
    }

    // If flight is enabled, let the particles evolve.
    if (s.fly) {
      for (auto &&p : s.particles) {
        auto &&xy = p.xy;
        xy += p.v * dt;
      }
      s.group();
    }
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
