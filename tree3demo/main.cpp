#include <algorithm>
#include <complex>
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
    // Compute groups/nodes.
    uint64_t constexpr MASK = 0xffff'ffff'ffff'0000;
    auto z_masked = [MASK](auto xy) {
      auto z = dyn::fixedmorton32<512>(xy);
      if (z.has_value())
        return std::optional<uint64_t>{z.value() & MASK};
      else
        return std::optional<uint64_t>{};
    };
    auto with_node = [&s](auto node) { s.nodes.push_back(node); };
    dyn::tree32::group<Physical>(s.particles.begin(), s.particles.end(),
                                 z_masked, with_node);
    for (auto &&n : s.nodes) {
      std::complex<float> center;
      float radius{}, count{1.0f};
      for (auto &&p : n)
        center += (p - center) / count++;
      for (auto &&p : n)
        radius = std::max(radius, std::abs(p - center));
      n.extra.center = center;
      n.extra.radius = radius;
    }
    return s;
  }

private:
  State() = default;
};

int do_main() {
  InitWindow(600, 600, "A");
  SetTargetFPS(15);
  auto s = State::fresh();
  Camera2D cam{};
  {
    auto w = float(GetScreenWidth()), h = float(GetScreenHeight()),
         z = std::min(w, h);
    cam.offset = {w * 0.5f, h * 0.5f};
    cam.zoom = 0.25f * z;
  }

  while (!WindowShouldClose()) {
    if (IsKeyPressed(KEY_R)) {
      s = State::fresh();
      BeginDrawing();
      EndDrawing();
      continue;
    }

    BeginDrawing();
    ClearBackground(BLACK);
    BeginMode2D(cam);
    {
      // Mouse.
      auto [mx, my] = GetScreenToWorld2D(GetMousePosition(), cam);
      auto mouse = std::complex{mx, my};

      // 2 px.
      auto radius = 2.0f / cam.zoom;
      for (auto &&c : s.nodes) {
        for (auto &&p : c) {
          DrawCircleV({p.real(), p.imag()}, radius, WHITE);
        }
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
