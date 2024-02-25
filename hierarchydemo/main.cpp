// Demonstrate Barnes-Hut approximation.

#include <algorithm>
#include <barnes_hut.h>
#include <cmath>
#include <complex>
#include <cstdint>
#include <numbers>
#include <optional>
#include <random>
#include <raylib.h>
#include <utility>
#include <vector>

#include "user.h"

/// Particle with location and velocity.
struct Particle {
  std::complex<float> xy, v;
  std::optional<uint64_t> morton_code{};

  Particle() = default;
  Particle(float x, float y, float vx, float vy) : xy{x, y}, v{vx, vy} {}

  void update_morton() { morton_code = dyn::bh32::morton(xy); }

  [[nodiscard]] std::optional<uint64_t> morton() const {
    return dyn::bh32::morton(xy);
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

  Physicals &operator+=(Physicals const &p) {
    if (this == &p)
      return count *= 2.0f, *this;
    typedef std::complex<double> C;
    typedef std::complex<float> F;
    auto sum = double(count) + double(p.count);
    auto proportion0 = double(count) / sum;
    auto proportion1 = double(p.count) / sum;
    count += p.count;
    xy = F{proportion0 * C{xy} + proportion1 * C{p.xy}};
    radius = std::max(radius, p.radius + std::abs(p.xy - xy));
    return *this;
  }
};

/// Store particles
struct State : public std::vector<Particle> {
  State(int N = 50'000) {
    std::mt19937 r(std::random_device{}());
    std::normal_distribution<float> z;
    for (auto i = 0; i < N; i++)
      emplace_back(z(r), z(r), z(r), z(r));
    sort();
  }

  void fly(float dt) {
    for (auto &&p : *this)
      p.xy += dt * p.v;
    sort();
  }

  void sort() {
    // Update the Morton codes.
    for (auto &&p : *this)
      p.update_morton();
    // Most particles stay where they used to be (if being called again).
    std::ranges::stable_sort(begin(), end(), {},
                             [](auto &&p) { return p.morton(); });
  }
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
        auto gxy = g.data.xy;
        auto gr = g.data.radius;
        auto dist = std::abs(gxy - w_mouse);
        // Eliminate the group (g) and all its descendant particles if the given
        // point (w_mouse) is far away from the boundary of the group's circle.
        if (max_view_distance < dist - gr)
          return ERASE;
        auto square = [](auto x) { return x * x; };
        auto dim = 1.0f - square(dist / max_view_distance); // [0, 1] closed.
        if (g.single())
          // This group (g) wraps a single particle.
          return user.dot(gxy, Fade(WHITE, dim)), ERASE;
        // Test the distance and the (approximate) viewing angle.
        if (dist < gr)
          // This group (g)'s circle contains the given point (w_mouse).
          // Higher level of detail required.
          return KEEP;
        // Construct a radius perpendicular to the line of sight from the
        // given point (w_mouse) from the center of the group's circle, and
        // then measure the angle between the ray from w_mouse to the radial
        // endpoint and the ray of the line of sight. This is an
        // underapproximation (but a good one) of one-half of the true view
        // angle.
        if (auto tan = gr / dist; tan_angle_threshold < tan)
          // View angle too wide; higher detail required.
          return KEEP;
        // View angle is small enough. Treat g as a point particle. Draw the
        // circle that represents g for visualization.
        DrawCircleLinesV({gxy.real(), gxy.imag()}, gr, Fade(YELLOW, dim));
        return ERASE;
      };

      namespace bh = dyn::bh32;

      // Require Z-sorted particles in state.
      bh::View<State const &> view{state};
      // Masked Morton (Z) code.
      auto morton = [](auto &&p, uint64_t m) -> std::optional<uint64_t> {
        if (auto w = p.morton(); w.has_value())
          return w.value() & m;
        else
          return {};
      };
      // Compute the "levels."
      auto levels = bh::levels<Physicals>(view, morton);
      bh::run(levels, process);
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
