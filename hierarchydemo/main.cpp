// Demonstrate Barnes-Hut approximation.

#include <barnes_hut.h>
#include <bit>
#include <cmath>
#include <complex>
#include <cstdio>
#include <list>
#include <morton.h>
#include <numbers>
#include <optional>
#include <random>
#include <raylib.h>
#include <vector>

/// Particle with location and velocity.
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

/// Used for signaling in the Barnes-Hut iteration.
enum { KEEP = 0, ERASE = 1 };

/// Barnes-Hut "group" of particles
template <typename I> struct Group {
  I first, last;
  std::complex<float> xy;
  float radius{};
  Group() = default;
  Group(I first, I last) : first(first), last(last) {}
  [[nodiscard]] I begin() { return first; }
  [[nodiscard]] I end() { return last; }
  [[nodiscard]] I begin() const { return first; }
  [[nodiscard]] I end() const { return last; }
  [[nodiscard]] bool single() const {
    if (first == last)
      return false;
    auto a{first};
    return ++a == last;
  }
};

/// Store particles
class State {
  /// A list of particles.
  std::vector<Particle> particles;
  friend class View;

public:
  State(int N = 50'000) {
    std::mt19937 r(std::random_device{}());
    std::normal_distribution<float> z;
    for (auto i = 0; i < N; i++)
      particles.emplace_back(z(r), z(r), z(r), z(r));
    sort();
  }

  void fly(float dt) {
    for (auto &&p : particles)
      p.xy += dt * p.v;
    sort();
  }

  void sort() {
    auto constexpr morton = [](Particle const &p) {
      return Particle::morton(p);
    };
    std::sort(particles.begin(), particles.end(),
              [&morton](auto &&a, auto &&b) { return morton(a) < morton(b); });
  }

  [[nodiscard]] size_t size() const { return particles.size(); }
};

/// Build groups, do Barnes-Hut stuff
class View {
  static uint64_t constexpr DEFAULT_MASK = 0xffff'ffff'0000'0000;
  State &s;
  uint64_t mask;

public:
  View(State &s) : s(s), mask(DEFAULT_MASK) {}

  using Groups = std::list<Group<decltype(s.particles.cbegin())>>;
  Groups groups(Groups &&prior = {}) {
    if (mask == uint64_t(-1) || s.particles.empty())
      return {};
    Groups r;
    auto m = mask;
    auto z = [m](auto &&p) { return Particle::morton(p, m); };
    auto grp = [&r](auto f, auto l) { r.push_back({f, l}); };
    if (prior.empty())
      // First time? Compute everything.
      dyn::bh32::group(s.particles.cbegin(), s.particles.cend(), z, grp);
    else
      // Subdivide group (g) knowing that no particle can jump through adjacent
      // groups by construction (sorted by Z-code, particles belong to disjoint
      // squares).
      for (auto &&g : prior)
        dyn::bh32::group(g.begin(), g.end(), z, grp);
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
  void reset_mask() { mask = 0xffff'ffff'0000'0000; }
  void shift_right() {
    auto m = std::bit_cast<int64_t>(mask);
    mask = std::bit_cast<uint64_t>(m >> 2);
  }
};

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

  void dot(auto &&p, auto color) {
    auto radius = 2.0f / cam.zoom;
    DrawCircleV({p.real(), p.imag()}, radius, color);
  }
};

static int do_main() {
  // Prepare the particles, and then Z-sort them.
  State state;

  SetConfigFlags(FLAG_WINDOW_RESIZABLE);
  InitWindow(600, 600, "X");
  SetTargetFPS(60);

  // (Auxiliary object for GUI).
  User user;

  // Performance is highly sensitive to value:
  auto const tan_angle_threshold =
      std::tan(8.0f * std::numbers::pi_v<float> / 180.0f);

  // Maximum viewing distance in world length units (or +infinity)
  auto const max_view_distance = 4.0f;

  while (!WindowShouldClose()) {
    BeginDrawing();
    ClearBackground(BLACK);
    BeginMode2D(user.cam);
    {
      user.pan(), user.zoom();
      if (user.adjust_fly(), user.flag.fly)
        state.fly(1.0f / 60.0f); // [1 / Hz] = [seconds].

      // Get an example "starter" point in world coordinates (w_mouse).
      auto [mx, my] = GetMousePosition();
      auto v_mouse = GetScreenToWorld2D({mx, my}, user.cam);
      auto w_mouse = std::complex{v_mouse.x, v_mouse.y};

      // Apply BFS (breadth-first search) to the implicit quadtree.

      // For each "group" (bunch of particles considered to have the same
      // level of detail), do whatever is desired, and then decide whether the
      // particles in the group should be considered at the next round (KEEP),
      // or if all the particles could be discarded (REMOVE).
      auto process = [&user, w_mouse, tan_angle_threshold,
                      max_view_distance](auto &&g) {
        auto dist = std::abs(g.xy - w_mouse);
        // Eliminate the group (g) and all its descendant particles if the given
        // point (w_mouse) is far away from the boundary of the group's circle.
        if (max_view_distance < dist - g.radius)
          return ERASE;
        auto square = [](auto x) { return x * x; };
        auto dim = 1.0f - square(dist / max_view_distance); // [0, 1] closed.
        if (g.single())
          // This group (g) wraps a single particle.
          return user.dot(g.xy, Fade(WHITE, dim)), ERASE;
        // Test the distance and the (approximate) viewing angle.
        if (dist < g.radius)
          // This group (g)'s circle contains the given point (w_mouse).
          // Higher level of detail required.
          return KEEP;
        // Construct a radius perpendicular to the line of sight from the
        // given point (w_mouse) from the center of the group's circle, and
        // then measure the angle between the ray from w_mouse to the radial
        // endpoint and the ray of the line of sight. This is an
        // underapproximation (but a good one) of one-half of the true view
        // angle.
        if (auto tan = g.radius / dist; tan_angle_threshold < tan)
          // View angle too wide; higher detail required.
          return KEEP;
        // View angle is small enough. Treat g as a point particle. Draw the
        // circle that represents g for visualization.
        DrawCircleLinesV({g.xy.real(), g.xy.imag()}, g.radius,
                         Fade(YELLOW, dim));
        return ERASE;
      };

      // Require Z-sorted particles in state.
      View view{state};
      for (decltype(view.groups()) groups;; view.shift_right()) {
        // [1] Construct quadtree nodes at given depth (the mask; stored in s).
        // If no groups exist (either because the depth limit has been reached
        // or because there are no more particles), terminate the loop.
        groups = view.groups(std::move(groups));
        // [2] Apply a linear scan to all the groups found at the depth. Process
        // the group. Erase the particles that must be erased. In the next
        // round, with an additional level of detail, there will be hopefully
        // fewer particles.
        std::erase_if(groups, process);
        if (groups.empty())
          break;
      }
      view.reset_mask();
    }
    EndMode2D();

    User::hud(int(state.size()));

    if (IsKeyPressed(KEY_R))
      state = {}, user = {};

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
