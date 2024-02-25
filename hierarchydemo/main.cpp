// Demonstrate Barnes-Hut approximation.

#include <algorithm>
#include <barnes_hut.h>
#include <bit>
#include <cmath>
#include <complex>
#include <numbers>
#include <optional>
#include <random>
#include <raylib.h>
#include <vector>

#include "user.h"

/// Particle with location and velocity.
struct Particle {
  std::complex<float> xy, v;
  std::optional<uint64_t> morton_code{};

  Particle() = default;
  Particle(float x, float y, float vx, float vy) : xy{x, y}, v{vx, vy} {}

  void update_morton() { morton_code = morton(xy); }

  [[nodiscard]] std::optional<uint64_t> morton_masked(uint64_t mask) const {
    if (morton_code.has_value())
      return morton_code.value() & mask;
    else
      return {};
  }
  static std::optional<uint64_t> morton(std::complex<float> p) {
    return dyn::bh32::morton(p);
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
};

/// Used for signaling in the Barnes-Hut iteration.
enum { KEEP = 0, ERASE = 1 };

struct Physicals {
  std::complex<float> xy;
  float radius{}, count{};

  Physicals() = default;
  template <class I> Physicals(I first, I last) {
    auto f = first;
    while (f != last)
      xy += (f++->xy - xy) / ++count;
    f = first;
    while (f != last)
      radius = std::max(radius, std::abs(f++->xy - xy));
  }

  template <class I> void merge(Physicals &p, I first, I last) {

  }
};

/// Store particles
class State {
  /// A list of particles.
  std::vector<Particle> particles;

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
    // Update the Morton codes.
    for (auto &&p : particles)
      p.update_morton();
    // At least on MSVC (as of writing), stable_sort is seen to be faster than
    // sort. I have no idea why. It just is.
    std::ranges::stable_sort(particles.begin(), particles.end(), {},
                             [](auto &&p) { return p.morton_code; });
  }

  [[nodiscard]] size_t size() const { return particles.size(); }
};

static int do_main() {
  SetConfigFlags(FLAG_WINDOW_RESIZABLE);
  InitWindow(600, 600, "X");
  SetTargetFPS(60);

  // Prepare the particles, and then Z-sort them.
  State state;

  // (Auxiliary object for GUI).
  User user;

  // Performance is highly sensitive to value (tangent of half of viewing angle)
  //  Smaller angle: bad for performance, ostensibly more "accurate"
  //  Larger angle: good for performance, less "accurate"
  // *accuracy: depends on application.
  auto const tan_angle_threshold =
      std::tan(5.0f * std::numbers::pi_v<float> / 180.0f);

  // Maximum viewing distance in world length units (or +infinity)
  auto const max_view_distance = 10.0f;

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
#define U [[maybe_unused]]
int WinMain(U void *_0, U void *_1, U void *_2, U int _3) { return do_main(); }
#else
int main() { return do_main(); }
#endif
