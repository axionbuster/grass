#include <barnes_hut.h>
#include <bit>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <list>
#include <morton.h>
#include <optional>
#include <random>
#include <raylib.h>
#include <vector>

struct Particle {
  std::complex<float> xy, v;

  Particle() = default;
  Particle(float x, float y, float vx, float vy) : xy{x, y}, v{vx, vy} {}

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

enum class Test { KEEP, REMOVE };

template <typename I> struct Group {
  I first, last;
  std::complex<float> xy;
  float radius{};
  Group() = default;
  Group(I first, I last) : first(first), last(last) {}
  I begin() { return first; }
  I end() { return last; }
  I begin() const { return first; }
  I end() const { return last; }
};

struct State {
  std::list<Particle> particles;
  uint64_t mask = 0xc000'0000'0000'0000;
  [[maybe_unused]] void shift_left() {
    mask = (mask <<= 2) ? mask : 0xc000'0000'0000'0000;
  }
  void shift_right() {
    auto m = std::bit_cast<int64_t>(mask);
    mask = std::bit_cast<uint64_t>(m >> 2);
  }

  void sort() {
    auto morton = [](Particle const &p) { return Particle::morton(p); };
    particles.sort(
        [&morton](auto &&a, auto &&b) { return morton(a) < morton(b); });
  }

  static State fresh(int N = 10'000) {
    decltype(particles) ps;
    std::mt19937 r(1234);
    std::normal_distribution<float> z;
    for (auto i = 0; i < N; i++)
      ps.emplace_back(z(r), z(r), z(r), z(r));
    State s{ps};
    s.sort();
    return s;
  }

  std::vector<Group<std::list<Particle>::iterator>> groups() {
    if (mask == uint64_t(-1) || particles.empty())
      return {};

    decltype(groups()) r;
    auto m = mask;
    auto z = [m](auto &&p) { return Particle::morton(p, m); };
    auto grp = [&r](auto f, auto l) { r.push_back({f, l}); };
    dyn::bh32::group(particles.begin(), particles.end(), z, grp);

    // Position, radius.
    for (auto &&g : r) {
      // Welford's online mean algorithm: compute mean xy for all particles
      // (p) in the group (g).
      float n{};
      for (auto &&p : g)
        g.xy += (p.xy - g.xy) / ++n;
      // Compute min radius to contain all particles (p).
      for (auto &&p : g)
        g.radius = std::max(g.radius, std::abs(p.xy - g.xy));
    }
    return r;
  }

  void fly(float dt) {
    for (auto &&p : particles)
      p.xy += dt * p.v;
    sort();
  }

private:
  State() = default;
  State(decltype(particles) ps) : particles{std::move(ps)} {}
};

struct User {
  Camera2D cam{};
  float zoom0{};
  struct {
    bool fly : 1 {};
  } flag{};

  User() {
    auto w = float(GetScreenWidth()), h = float(GetScreenHeight()),
         z = std::min(w, h);
    cam.offset = {w * 0.5f, h * 0.5f};
    cam.zoom = 0.25f * z;
    zoom0 = cam.zoom;
  }

  void pan() {
    if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) {
      auto u = GetMouseDelta();
      auto v = std::complex{u.x, u.y};
      auto w = std::complex{cam.target.x, cam.target.y};
      auto x = w - v / cam.zoom;
      cam.target = {x.real(), x.imag()};
    }
  }

  void zoom() {
    if (auto wheel = GetMouseWheelMove()) {
      auto u = GetMousePosition();
      auto v = GetScreenToWorld2D(u, cam);
      cam.offset = u;
      cam.target = v;
      auto constexpr ZOOM_INCR = 5.0f;
      cam.zoom += wheel * ZOOM_INCR;
      cam.zoom = std::max(cam.zoom, 0.25f * zoom0);
      cam.zoom = std::min(cam.zoom, 10.0f * zoom0);
    }
  }

  void hud(int n_particles) {
    DrawFPS(16, 16);
    char msg[256]{};
    snprintf(msg, sizeof msg, "%d particles", n_particles);
    DrawText(msg, 16, 40, 20, WHITE);
  }

  void adjust_fly() {
    if (IsKeyPressed(KEY_SPACE))
      flag.fly = !flag.fly;
  }
};

static void draw_particle(auto &&p, auto color, auto zoom) {
  auto radius = 2.0f / zoom;
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
    u.hud(int(state.particles.size()));
    BeginMode2D(u.cam);
    {
      u.pan(), u.zoom();
      if (u.adjust_fly(), u.flag.fly)
        // dt (Hz)
        state.fly(1.0f / 60.0f);

      // I'll explain everything soon when I get the time to work on this again.
      auto s = state;
      auto [mx, my] = GetMousePosition();
      auto v_mouse = GetScreenToWorld2D({mx, my}, u.cam);
      auto w_mouse = std::complex{v_mouse.x, v_mouse.y};
      auto tan_angle_threshold = 0.176327f; // tan (10 deg).
      for (decltype(s.groups()) groups; !(groups = s.groups()).empty();
           s.shift_right()) {
        auto process = [&u, w_mouse, tan_angle_threshold](auto &&g) {
          if (!g.radius) {
            draw_particle(g.xy, WHITE, u.cam.zoom);
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
        for (auto &&g : groups)
          if (process(g) == Test::REMOVE)
            s.particles.erase(g.begin(), g.end());
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
