#ifndef GRASS_USER_H
#define GRASS_USER_H

#include <complex>
#include <iosfwd>
#include <optional>
#include <raylib.h>
#include <utility>

/// User interface
struct User {
  /// Current camera (Raylib).
  Camera2D cam{};

  /// Control over features.
  struct {
    bool fly : 1 {true};
    bool demo : 1 {true};
    // If set, skip spawning in this frame when user asks to spawn a particle.
    bool spawned_last_frame : 1 {};
    // Target FPS.
    unsigned short target_fps{90};
    // Last time the simulation began or reset, seconds.
    double last_sec{GetTime()};
    [[nodiscard]] float target_dt() const { return 1.0f / float(target_fps); }
  } control;

  /// Show options (for the HUD).
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

  /// Original zoom level.
  float zoom0{};

  /// Construct a new instance.
  User() {
    auto w = float(GetScreenWidth()), h = float(GetScreenHeight()),
         z = std::min(w, h);
    cam.offset = {w * 0.5f, h * 0.5f};
    cam.zoom = 0.125f * z;
    zoom0 = cam.zoom;
  }

  /// If R is pressed, the user wants to reset the simulation.
  bool wants_reset() const { return IsKeyPressed(KEY_R); }

  /// If T is pressed, show a different option.
  void rotate_debug_opts() {
    if (IsKeyPressed(KEY_T))
      show.next();
  }

  /// Handle pan on input.
  void pan() {
    if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) {
      auto u = GetMouseDelta();
      auto v = std::complex{u.x, u.y};
      auto w = std::complex{cam.target.x, cam.target.y};
      auto x = w - v / cam.zoom;
      cam.target = {x.real(), x.imag()};
    }
  }

  /// Handle zoom on input.
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

  /// Measure the time since last reset.
  double elapsed_sec() const { return GetTime() - control.last_sec; }

  /// Test whether the user wants to spawn a particle right now. If so, return
  /// where (in world coordinates).
  std::optional<std::complex<float>> wants_spawn_particle() const {
    if (IsMouseButtonDown(MOUSE_BUTTON_LEFT))
      return mouse();
    else
      return {};
  }

  /// Return the mouse position in world coordinates.
  std::complex<float> mouse() const {
    auto a = GetScreenToWorld2D(GetMousePosition(), cam);
    return {a.x, a.y};
  }

  /// Write text.
  void hud(auto n_particles) const {
    // The standard library understands how to format a complex number, but,
    // understandably, knows nothing about Raylib's custom vector types.
    auto constexpr v2c = [](Vector2 v) {
      return std::complex<float>{v.x, v.y};
    };
    std::stringstream buf;
    if (control.demo)
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

  /// If SPACE is pressed, toggle flight.
  void adjust_fly() {
    if (IsKeyPressed(KEY_SPACE))
      control.fly = !control.fly;
  }

  /// Draw a particle.
  void particle(auto p, Color color = WHITE) const {
    DrawCircleV({p.real(), p.imag()}, p.radius, color);
  }

  /// A rectangle with the less-less and greater-greater coordinates.
  struct ComplexRectangle {
    /// Less-less (ll) and greater-greater (gg) coordinates.
    std::complex<float> ll, gg;
  };

  /// Get the rectangle (in world coordinates) that represents the window.
  [[nodiscard]] ComplexRectangle window() const {
    auto w = float(GetScreenWidth()), h = float(GetScreenHeight());
    auto a = GetScreenToWorld2D({}, cam), b = GetScreenToWorld2D({w, h}, cam);
    return {{a.x, a.y}, {b.x, b.y}};
  }
};

#endif /* GRASS_USER_H */
