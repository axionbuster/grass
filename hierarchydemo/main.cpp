// Demonstrate Barnes-Hut approximation.

#include <algorithm>
#include <barnes_hut.h>
#include <bit>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <list>
#include <morton.h>
#include <numbers>
#include <optional>
#include <random>
#include <raylib.h>
#include <type_traits>
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
  /// Test: Are the iterators equal (discarding xy and radius)?
  [[nodiscard]] bool operator==(Group const &g) const {
    if (this == &g)
      return true;
    // Ignore: xy, radius.
    return first == g.first && last == g.last;
  }
};

/// Store particles
class State {
  /// A list of particles.
  std::vector<Particle> particles;
  friend class View;

public:
  State(int N = 50'000) {
    std::mt19937 r(1234);
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

/// Build groups, do Barnes-Hut stuff
class View {
  State const &s;

  // Start fully general (from the lowest level of detail).
  //  (0xffff....ffff is the highest level of detail).
  uint64_t mask = 0xc000'0000'0000'0000;

public:
  /// Take an immutable view to a table of particles.
  /// Start with the lowest level of detail
  /// (call `refine` to shift gear to a higher level of detail).
  View(State const &s) : s(s) {}

  using Groups = std::vector<Group<decltype(s.particles.cbegin())>>;

  /// Compute the groups at the given level of detail, reusing work from an
  /// earlier call by the same instance of View for the same instance of State.
  /// (If not given, as usual for a fresh start, then compute everything).
  Groups groups(Groups const &prior = {}) const {
    if (!mask || s.particles.empty())
      // !mask <-> halt. Sentinel used by `refine()`.
      return {};
    Groups novel;
    auto m = mask;
    auto z = [m](auto &&p) { return p.morton_masked(m); };
    auto grp = [&novel](auto f, auto l) { novel.push_back({f, l}); };
    if (prior.empty())
      // First time? Compute everything.
      dyn::bh32::group(s.particles.cbegin(), s.particles.cend(), z, grp);
    else
      // Subdivide group (g) knowing that no particle can jump through adjacent
      // groups by construction (sorted by Z-code, particles belong to disjoint
      // squares).
      for (auto &&g : prior)
        dyn::bh32::group(g.begin(), g.end(), z, grp);
    // Test: same nodes pointing to the same particles? (ignoring physical
    // summary).
    if (prior == novel)
      // Since loop below is hot, if point to same particles, reuse the work.
      // Return `prior` which has the physical summary calculated.
      // Don't return `novel`: it has zeroes instead of actual physical data.
      return prior;
    // Compute physical summary.
    // Position, radius.
    for (auto &&g : novel) {
      // Welford's online mean algorithm: compute mean xy for all particles
      // (p) in the group (g).
      float n{};
      for (auto &&p : g)
        g.xy += (p.xy - g.xy) / ++n;
      // Compute min radius to contain all particles (p).
      for (auto &&p : g)
        g.radius = std::max(g.radius, std::abs(p.xy - g.xy));
    }
    return novel;
  }

  /// Shift to the next (higher) level of detail.
  void refine() {
    if (mask == 0xffff'ffff'ffff'ffff)
      // Signal `groups` to produce nothing and halt.
      mask = 0;
    // Apply arithmetic (sign-bit-extending) bit shift by two bits.
    // (Z-codes (aka Morton codes) use two bits to designate each quadrant.)
    auto m = std::bit_cast<int64_t>(mask);
    mask = std::bit_cast<uint64_t>(m >> 2);
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

  // (Re-use memory when computing groups).
  std::vector<View::Groups> levels;

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

      // For a realistic scenario, generate all the levels at once.
      levels.clear();
      levels.emplace_back(view.groups());
      do {
        // As expected, presents the greatest bottleneck.
        levels.emplace_back(view.groups(levels.back()));
        view.refine();
      } while (!levels.back().empty());
      levels.pop_back();

      // Remove adjacent duplicate layers.
      // (Ordered by inclusion; elements closer to begin() are either supersets
      // of or disjoint to the elements closer to end()).
      auto last = std::unique(levels.begin(), levels.end());

      // Traverse on copies. (Required to do so because it's a prototype to the
      // real one, which requires the tree to be built once per frame and then
      // traversed from top level to bottom level by copy.) Work in adjacent
      // pairs of layers from the first to last. In each pair of layers, only
      // accept the groups (elements of the layers) that are subsets of the
      // other layer's groups, and trim the rest.
      if (levels.begin() != last) {
        auto prior{*levels.begin()};
        auto begin = levels.begin();

        // Re-use allocated memory.
        decltype(prior) copy;
        while (++begin != last) {
          copy.clear();
          auto const &novel = *begin;
          auto pg = prior.begin();
          auto ng = novel.begin();
          if (ng == novel.end() || pg == prior.end())
            break;
          // Skip.
          while (ng->begin() != pg->begin())
            ++ng;
          // Cull.
          while (pg != prior.end()) {
            if (ng->begin() == pg->begin()) {
              while (ng->end() != pg->end())
                copy.push_back(*ng++);
              copy.push_back(*ng++), ++pg;
            } else
              while (ng->begin() != pg->begin())
                ++ng;
          }
          // Process.
          std::erase_if(copy, process);
          // Break if no more.
          if (copy.empty())
            break;
          prior = copy;
        }
      }
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
