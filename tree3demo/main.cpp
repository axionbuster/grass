#include <algorithm>
#include <bit>
#include <cinttypes>
#include <complex>
#include <cstdio>
#include <morton.h>
#include <optional>
#include <random>
#include <raylib.h>
#include <tree.h>
#include <vector>

int constexpr N = 1000;

struct Physical {
  std::complex<float> center;
  float radius{};
};

struct State {
  std::vector<std::complex<float>> particles;
  std::vector<dyn::tree32::Node<Physical, decltype(particles.begin())>> nodes;
  int64_t mask = 0xffff'ffff'ffff'0000;
  static State fresh() {
    State s;
    // Make random particles.
    std::mt19937 rng(std::random_device{}());
    std::normal_distribution<float> dist;
    s.particles.reserve(N);
    for (int i = 0; i < N; i++)
      s.particles.emplace_back(dist(rng), dist(rng));
    // Sort by morton.
    std::ranges::sort(s.particles.begin(), s.particles.end(), {},
                      dyn::fixedmorton32<512>);
    s.group();
    return s;
  }

  void group() {
    nodes = {};
    // Compute groups/nodes.
    auto m = std::bit_cast<uint64_t>(mask);
    auto z_masked = [m](auto xy) {
      auto z = dyn::fixedmorton32<512>(xy);
      if (z.has_value())
        return std::optional<uint64_t>{z.value() & m};
      else
        return std::optional<uint64_t>{};
    };
    auto with_node = [this](auto node) { nodes.push_back(node); };
    dyn::tree32::group<Physical>(particles.begin(), particles.end(), z_masked,
                                 with_node);
    // Assign the extra data.
    for (auto &&n : nodes) {
      std::complex<float> center;
      float radius{}, count{};
      for (auto &&p : n)
        // Welford's online mean.
        center += (p - center) / ++count;
      for (auto &&p : n)
        radius = std::max(radius, std::abs(p - center));
      n.extra.center = center;
      n.extra.radius = radius;
    }
  }

  // Each quadrant requires two bits.

  void mask_left() {
    if (!(mask <<= 2)) {
      mask = 0xc000'0000'0000'0000;
    };
  }

  void mask_right() { mask >>= 2; }

private:
  State() = default;
};

int do_main() {
  SetConfigFlags(FLAG_WINDOW_RESIZABLE);
  InitWindow(600, 600, "A");
  SetTargetFPS(60);
  auto s = State::fresh();
  Camera2D cam{};
  float shortest_zoom;
  {
    // Do once at startup.
    auto w = float(GetScreenWidth()), h = float(GetScreenHeight()),
         z = std::min(w, h);
    cam.offset = {w * 0.5f, h * 0.5f};
    cam.zoom = 0.25f * z;
    shortest_zoom = cam.zoom;
  }

  while (!WindowShouldClose()) {
    if (IsKeyPressed(KEY_R)) {
      s = State::fresh();
      BeginDrawing();
      EndDrawing();
      continue;
    }

    if (IsKeyPressed(KEY_LEFT)) {
      s.mask_left();
      s.group();
    }

    if (IsKeyPressed(KEY_RIGHT)) {
      s.mask_right();
      s.group();
    }

    if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) {
      auto u = GetMouseDelta();
      auto v = std::complex{u.x, u.y};
      auto w = std::complex{cam.target.x, cam.target.y};
      auto x = w - v / cam.zoom;
      cam.target = {x.real(), x.imag()};
    }

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
    {
      auto mask = std::bit_cast<uint64_t>(s.mask);
      char text[128]{};
      std::snprintf(text, sizeof text, "Mask = 0x%" PRIx64, mask);
      DrawText(text, 16, 16, 20, YELLOW);
      std::snprintf(text, sizeof text, "(%d nodes/groups)",
                    int(s.nodes.size()));
      DrawText(text, 16, 40, 20, YELLOW);
      DrawText("<Left or right key to shift mask>", 16, 64, 20, YELLOW);
    }
    BeginMode2D(cam);
    {
      // Mouse.
      auto [mx, my] = GetScreenToWorld2D(GetMousePosition(), cam);
      auto mouse = std::complex{mx, my};

      // 2 px.
      auto radius = 2.0f / cam.zoom;
      for (auto &&c : s.nodes) {
        for (auto &&p : c)
          DrawCircleV({p.real(), p.imag()}, radius, WHITE);
        auto &&e = c.extra;
        if (e.radius) {
          auto center = Vector2{e.center.real(), e.center.imag()};
          auto distance = std::abs(mouse - e.center);
          if (distance < e.radius) {
            auto fade = 1.0f - distance / e.radius;
            auto color = Fade(WHITE, fade);
            DrawCircleV(center, e.radius, color);
          }
          DrawCircleLinesV(center, e.radius, WHITE);
        }
      }
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
