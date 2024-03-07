#include <cmath>
#include <complex>
#include <concepts>
#include <cstdio>
#include <raylib.h>
#include <sstream>

/// User interface
struct User {
  Camera2D cam{};
  float zoom0{};
  struct {
    bool fly : 1 {true};
    // If set, skip spawning in this frame when user asks to spawn a particle.
    bool spawned_last_frame : 1 {};
    // Target FPS.
    unsigned short target_fps{90};
    // Last time the simulation began or reset, seconds.
    double last_sec{GetTime()};
    explicit operator bool() const { return fly; }
    [[nodiscard]] float target_dt() const { return 1.0f / float(target_fps); }
  } control;
  struct {
    bool fps : 1 {};
    bool n_particles : 1 {};
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

  User() {
    auto w = float(GetScreenWidth()), h = float(GetScreenHeight()),
         z = std::min(w, h);
    cam.offset = {w * 0.5f, h * 0.5f};
    cam.zoom = 0.125f * z;
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

  void hud(int n_particles) const {
    // The standard library understands how to format a complex number, but,
    // understandably, knows nothing about Raylib's custom vector types.
    auto constexpr v2c = [](Vector2 v) {
      return std::complex<float>{v.x, v.y};
    };
    std::stringstream buf;
    if (!control)
      buf << "(Demo; click anywhere to add particles.)\n";
    if (!show.any())
      buf << "R to reset; T for debug";
    if (show.fps)
      buf << "FPS: " << GetFPS() << '\n';
    if (show.n_particles)
      buf << "N: " << n_particles << '\n';
    if (show.cam)
      buf << "Zoom: " << cam.zoom << "\nTarget: " << v2c(cam.target)
          << "\nOffset: " << v2c(cam.offset) << '\n';
    auto str = buf.str();
    DrawText(str.c_str(), 16, 16, 20, LIGHTGRAY);
  }

  void adjust_fly() {
    if (IsKeyPressed(KEY_SPACE))
      control.fly = !control.fly;
  }

  void particle(auto p, auto color) const {
    DrawCircleV({p.real(), p.imag()}, p.radius, color);
  }

  struct ComplexRectangle {
    std::complex<float> ll, gg;
  };

  [[nodiscard]] ComplexRectangle window_world() const {
    auto uv = std::complex{cam.offset.x, cam.offset.y} / cam.zoom;
    auto target = std::complex{cam.target.x, cam.target.y};
    return {-uv + target, uv + target};
  }
};
