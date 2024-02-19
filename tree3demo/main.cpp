#include <complex>
#include <morton.h>
#include <optional>
#include <random>
#include <ranges>
#include <raylib.h>
#include <tree.h>
#include <vector>

int constexpr N = 1000;

struct State {
  std::vector<std::complex<float>> particles;
  std::vector<
      dyn::tree32::Node<std::ranges::subrange<decltype(particles.begin())>>>
      wip;
  State() {
    std::mt19937 rng(std::random_device{}());
    std::normal_distribution<float> dist;
    particles.reserve(N);
    for (int i = 0; i < N; i++)
      particles.emplace_back(dist(rng), dist(rng));
  }
};

template <std::ranges::range R> struct MyRange {
  R range;
  MyRange(auto begin, auto end) : range{begin, end} {}
};

int do_main() {
  InitWindow(600, 600, "A");
  SetTargetFPS(60);
  State s;
  Camera2D cam{};
  {
    auto w = float(GetScreenWidth()), h = float(GetScreenHeight()),
         z = std::min(w, h);
    cam.offset = {w * 0.5f, h * 0.5f};
    cam.zoom = 0.25f * z;
  }
  while (!WindowShouldClose()) {
    {
      // Compute morton.
      uint64_t constexpr MASK = 0xffff'ffff'ffff'0000;
      auto z_masked = [MASK](auto xy) {
        auto z = dyn::fixedmorton32<512>(xy);
        if (z.has_value())
          return z.value() & MASK;
        else
          return std::optional<uint64_t>{};
      };
      //      auto with_node = [&s](auto );
    }

    BeginDrawing();
    ClearBackground(BLACK);
    BeginMode2D(cam);
    {
      // 2 px.
      auto radius = 2.0f / cam.zoom;
      auto &&pp = s.particles;
      for (auto &&p : pp)
        DrawCircleV({p.real(), p.imag()}, radius, WHITE);
      // Tree circles.
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
