#include <algorithm>
#include <barnes_hut.h>
#include <bit>
#include <circle.h>
#include <cmath>
#include <complex>
#include <cstdint>
#include <morton.h>
#include <optional>
#include <random>
#include <raylib.h>
#include <vector>

struct Particle {
  std::complex<float> xy;

  Particle() = default;
  Particle(float x, float y) : xy{x, y} {}

  static std::optional<uint64_t> morton(std::complex<float> p) {
    return dyn::fixedmorton32(p);
  }

  static std::optional<uint64_t> morton(Particle const &p) {
    return morton(p.xy);
  }

  static std::optional<uint64_t> morton(std::complex<float> p, uint64_t mask) {
    if (auto z = morton(p); z.has_value())
      return z.value() & mask;
    else
      return {};
  }

  static std::optional<uint64_t> morton(Particle const &p, uint64_t mask) {
    return morton(p.xy, mask);
  }
};

template <typename I> struct Group {
  I first, last;
  std::complex<float> xy;
  float radius{};
  Group() = default;
  Group(I first, I last) : first(first), last(last) {}
};

struct State {
  std::vector<Particle> particles;
  uint64_t mask = 0xffff'ffff'ffff'0000;
  void shift_left() { mask = (mask <<= 2) ? mask : 0xc000'0000'0000'0000; }
  void shift_right() {
    // Right shift on signed integers was undefined until C++20.
    // Since C++20, it's defined as an arithmetic (zero-extending) right shift.
    auto m = std::bit_cast<int64_t>(mask);
    mask = std::bit_cast<uint64_t>(m >> 2);
  }
  static State fresh() {
    decltype(particles) ps;
    auto constexpr N = 1'000;
    std::mt19937 r(1234);
    std::normal_distribution<float> z;
    for (auto i = 0; i < N; i++)
      ps.emplace_back(z(r), z(r));
    auto morton = [](Particle const &p) { return Particle::morton(p); };
    std::ranges::sort(ps.begin(), ps.end(), {}, morton);
    return {ps};
  }

  std::vector<Group<std::vector<Particle>::iterator>> groups() {
    if (mask + 1 == 0 || particles.empty())
      return {};

    decltype(groups()) r;
    auto m = mask;
    auto z = [m](auto &&p) { return Particle::morton(p, m); };
    auto g = [&r](auto f, auto l) { r.push_back({f, l}); };
    dyn::bh32::group(particles.begin(), particles.end(), z, g);
    return r;
  }

private:
  State() = default;
  State(decltype(particles) ps) : particles{std::move(ps)} {}
};

struct User {
  Camera2D cam;

  User() {
    auto w = float(GetScreenWidth()), h = float(GetScreenHeight()),
         z = std::min(w, h);
    cam.offset = {w * 0.5f, h * 0.5f};
    cam.zoom = 0.25f * z;
  }
};

static void draw_particle(auto &&p, auto color, auto &&cam) {
  auto radius = 2.0f / cam.zoom;
  DrawCircleV({p.real(), p.imag()}, radius, color);
}

static int do_main() {
  auto state = State::fresh();
  SetConfigFlags(FLAG_WINDOW_RESIZABLE);
  InitWindow(600, 600, "X");
  SetTargetFPS(60);
  User u;
  while (!WindowShouldClose()) {
    BeginDrawing();
    ClearBackground(BLACK);
    BeginMode2D(u.cam);
    {
      auto s = state;

      if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
        auto [mx, my] = GetMousePosition();
        auto v_mouse = GetScreenToWorld2D({mx, my}, u.cam);
        auto w_mouse = std::complex{v_mouse.x, v_mouse.y};
        auto tan_angle_threshold = 0.1763270f; // tan (10 deg).
        for (auto groups = s.groups(); !groups.empty();) {
          auto process = [&u, w_mouse, tan_angle_threshold](auto &&g) {
            using Test = dyn::bh32::Test;

            if (!g.radius) {
              draw_particle(g.xy, WHITE, u.cam);
              return Test::REMOVE;
            }

            auto disp = g.xy - w_mouse;
            auto dist = std::abs(disp);
            if (dist < g.radius)
              goto keep;
            if (auto tan = g.radius / dist; tan_angle_threshold < tan)
              goto keep;

            DrawCircleLinesV({g.xy.real(), g.xy.imag()}, g.radius, YELLOW);
            return Test::REMOVE;

          keep:
            DrawCircleLinesV({g.xy.real(), g.xy.imag()}, g.radius, GRAY);
            return Test::KEEP;
          };
          auto remove = [&s](auto &&g) { s.particles.erase(g.first, g.last); };
          dyn::bh32::run_level(groups.rbegin(), groups.rend(), process, remove);
          s.shift_right();
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
int WinMain(void *_0, void *_1, void *_2, int _3) { return do_main(); }
#else
int main() { return do_main(); }
#endif
