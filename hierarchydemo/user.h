#include <cmath>
#include <complex>
#include <raylib.h>

/// User interface
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

  static void hud(int n_particles) {
    // px (screen); offset (16) + [font size {20} + padding {4}] * line.
    auto y = [](int i) { return 16 + 24 * i; };
    auto i = 0;
    DrawFPS(16, y(i++));
    char msg[64]{};
    snprintf(msg, sizeof msg, "%d particles", n_particles);
    DrawText(msg, 16, y(i++), 20, WHITE);
    auto dt_ms = 1000.0f * GetFrameTime();
    snprintf(msg, sizeof msg, "time %.1f ms", dt_ms);
    DrawText(msg, 16, y(i++), 20, WHITE);
  }

  void adjust_fly() {
    if (IsKeyPressed(KEY_SPACE))
      flag.fly = !flag.fly;
  }

  void dot(auto &&p, auto color) const {
    auto radius = 2.0f / cam.zoom;
    DrawCircleV({p.real(), p.imag()}, radius, color);
  }
};
