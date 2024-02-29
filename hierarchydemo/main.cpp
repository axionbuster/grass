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

/// Particle with position and velocity.
struct Particle {
  /// Position and velocity
  std::complex<float> xy, v;

  /// Since the Morton (Z) code calculation has a huge overhead
  /// (thought to be mostly due to floating-point-to-integer conversion), it's
  /// good to precompute them.
  std::optional<uint64_t> morton_code{};

  Particle() = default;

  Particle(float x, float y, float vx, float vy) : xy{x, y}, v{vx, vy} {}

  /// Compute the Morton (Z) code and save it. The internal morton code may
  /// still have no value if a floating-point issue occurs (for instance,
  /// value too big).
  void update_morton() { morton_code = dyn::bh32::morton(xy); }

  /// Recall the last stored Morton (Z) code.
  [[nodiscard]] std::optional<uint64_t> morton() const { return morton_code; }
};

/// "Extra data" stored for a Barnes-Hut tree node. A circle.
struct Physicals {
  /// Center
  std::complex<float> xy;

  /// Radius and number of particles (integer).
  /// (Integers in floating points to avoid integer-float conversions).
  float radius{}, count{};

  // Three required member functions (by `dyn::bh32::View::layer`):
  //  1. No-argument constructor.
  //  2. Constructor given a range of particles.
  //  3. Merger function using (+=).

  Physicals() = default;

  /// Given a range of particles (with an `xy` field), compute the quantities.
  template <class I> Physicals(I first, I last) {
    auto f = first;
    while (f != last)
      // Welford's online average algorithm.
      xy += (f++->xy - xy) / ++count;
    f = first;
    while (f != last)
      radius = std::max(radius, std::abs(f++->xy - xy));
  }

  /// Merge p's information.
  Physicals &operator+=(Physicals const &p) {
    if (this == &p)
      // Cannot violate const contract.
      return count *= 2.0f, *this;
    // Compute the new average xy.
    auto sum = count + p.count;
    auto proportion0 = count / sum;
    auto proportion1 = p.count / sum;
    xy = proportion0 * xy + proportion1 * p.xy;
    // The rest.
    count += p.count;
    radius = std::max(radius, p.radius + std::abs(p.xy - xy));
    return *this;
  }

  [[nodiscard]] bool single() const { return count == 1.0f; }
};

/// Store particles
struct State : public std::vector<Particle> {
  /// Construct a few particles.
  State(int N = 5) {
    std::mt19937 r(std::random_device{}());
    std::normal_distribution<float> z;
    for (auto i = 0; i < N; i++)
      emplace_back(z(r), z(r), z(r), z(r));
    sort();
  }

  /// Let the particles fly linearly for a while.
  /// @param dt Time in unit time [T].
  void fly(float dt) {
    for (auto &&p : *this)
      p.xy += dt * p.v;
    sort();
  }

  /// If particle position have been manually edited, call this function to
  /// maintain the invariant that the particles are Z-ordered.
  void sort() {
    // Update the Morton codes.
    for (auto &&p : *this)
      p.update_morton();
    // Get the Morton (or "Z") codes.
    auto z = [](auto &&p) { return p.morton(); };
    // Most particles stay where they used to be (if being called again).
    std::ranges::stable_sort(begin(), end(), {}, z);
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

      // Used for signaling in the Barnes-Hut iteration.
      // Contract (boolean-convertible) is due to the free function
      // `dyn::bh32::run`.
      enum { KEEP = 0, ERASE = 1 };

      // For each "group" (bunch of particles considered to have the same
      // level of detail), do whatever is desired, and then decide whether the
      // particles in the group should be considered at the next round (KEEP),
      // or if all the particles could be discarded (REMOVE).
      auto process = [&user, w_mouse, tan_angle_threshold,
                      max_view_distance](auto &&g) {
        auto dist = std::abs(g.xy - w_mouse);
        if (max_view_distance < dist - g.radius)
          // Too far from the boundary of the group's circle.
          return ERASE;
        auto square = [](auto x) { return x * x; };
        auto dim = 1.0f - square(dist / max_view_distance); // [0, 1] closed.
        if (g.single())
          // This group (g) wraps a single particle.
          // Process and then forget.
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

      namespace bh = dyn::bh32;

      // Masked Morton (Z) code.
      auto morton = [](auto &&p, uint64_t m) -> std::optional<uint64_t> {
        if (auto w = p.morton(); w.has_value())
          return w.value() & m;
        else
          return {};
      };
      auto tree = bh::tree<Physicals>(state.begin(), state.end(), morton);
      bh::run(&*tree, process);
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
