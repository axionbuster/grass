#include <algorithm>
#include <array>
#include <deque>
#include <morton.h>
#include <random>
#include <raylib.h>
#include <vector>

std::array<float, 3> hsl2rgb(std::array<float, 3> hsl) {
  // From Wikipedia. Don't ask me what these variables stand for or how this
  // works, but you can correspond these letters and notation to the math.
  auto l = hsl[2], a = hsl[1] * std::min(l, 1.0f - l);
  auto f = [l, a, &hsl](auto n) {
    auto k = std::remainder(n + hsl[0] * 12.0f, 12.0f);
    auto m = std::min(k - 3.0f, std::min(9.0f - k, 1.0f));
    return l - a * std::max(-1.0f, m);
  };
  return {f(0.0f), f(8.0f), f(4.0f)};
}

/// @brief Do main loop. Return `true` if to be shown again.
bool show() {
  long constexpr N_PARTICLES = 5000, MAX_N_QUEUE = N_PARTICLES / 4;
  float constexpr RADIUS = 0.0078125f, SATURATION = 0.75f, LIGHTNESS = 0.66f;
  uint64_t constexpr MASK = 0xffff'ffff'ffff'0000;
  auto const MORTON = dyn::fixedmorton32<512>;

  std::random_device seed;
  std::mt19937 rng(seed());
  std::uniform_real_distribution<float> udist;

  // Generate particles and then sort them by Morton order.
  std::vector<std::complex<float>> pp;
  for (long i = 0; i < N_PARTICLES; i++) {
    auto point = [&rng, &udist]() { return 4.0f * udist(rng) - 2.0f; };
    pp.emplace_back(point(), point());
  }
  std::ranges::sort(pp.begin(), pp.end(), {}, MORTON);

  // Let center of window point to world origin.
  auto scw = float(GetScreenWidth()), sch = float(GetScreenHeight());
  auto scs = std::min(scw, sch) * 0.5f;
  Camera2D cam{};
  cam.zoom = scs * 0.5f;
  cam.offset = {scs, scs};

  // Let the drawing begin at the particle at index q.
  long q{};
  // Morton prefix (same color or different color?) See also: MASK.
  uint64_t prefix{};
  std::deque<Color> colors;
  while (!WindowShouldClose()) {
    BeginDrawing();
    // If `R` pressed, that means to reset.
    // Return while returning `true` (do this again).
    if (IsKeyPressed(KEY_R))
      // Return after `BeginDrawing` to avoid freeze.
      // (This also means to put `EndDrawing` here as well.)
      return EndDrawing(), true;
    BeginMode2D(cam);
    {
      ClearBackground(BLACK);
      // Color and then draw the particles.
      auto particle = pp[q];
      auto morton = MORTON(particle).value_or(prefix);
      // Compare prefix with previous particle to decide new color.
      auto pfx = morton & MASK;
      if (pfx == prefix) {
        // Same quadrant, old color.
        colors.push_front(colors.front());
        if (colors.size() > MAX_N_QUEUE)
          colors.pop_back();
      } else {
        // New quadrant, new color.
        auto hue = udist(rng);
        auto rgb = hsl2rgb({hue, SATURATION, LIGHTNESS});
        for (auto &&c : rgb)
          c *= 255.0f;
        typedef unsigned char Q;
        Color c{Q(rgb[0]), Q(rgb[1]), Q(rgb[2]), 255};
        colors.push_front(c);
        if (colors.size() > MAX_N_QUEUE)
          colors.pop_back();
      }
      prefix = pfx;
      // Draw particle at index q but also the trailing particles (n).
      for (long n = 0; n < MAX_N_QUEUE; n++) {
        // Compute actual index (j) of the particle (p) as a function of index q
        // and trail number n.
        auto i = q - n;
        auto j = i < 0 ? i + N_PARTICLES : i;
        auto p = pp[j];
        // Manipulate color (fade out).
        auto c = n < long(colors.size()) ? colors[n] : BLACK;
        auto ease_cubic = [](float a) { return 1.0f - a * a * a; };
        c = Fade(c, ease_cubic(float(n) / float(MAX_N_QUEUE)));
        // Yup. Draw it.
        DrawCircleV({p.real(), p.imag()}, RADIUS, c);
      }
      q = (q + 1) % N_PARTICLES;
      {
        auto code = dyn::fixedmorton32<512>(pp[0]).value() & MASK;
        for (;;) {
          auto proj = [MASK](auto c) -> std::optional<uint64_t> {
            if (auto d = dyn::fixedmorton32<512>(c); d.has_value()) {
              return d.value() & MASK;
            } else {
              return {};
            }
          };
          auto [first, last] =
              std::ranges::equal_range(pp.begin(), pp.end(), code, {}, proj);
          if (first == last)
            break;
          last--;
          auto ll = *first, gg = *last;
          auto wh = gg - ll;
          Rectangle rect{ll.real(), ll.imag(), wh.real(), wh.imag()};
          DrawRectangleLinesEx(rect, RADIUS, WHITE);
          last++;
          if (last == pp.end())
            break;
          if (auto c = dyn::fixedmorton32<512>(*last); c.has_value()) {
            code = c.value() & MASK;
          } else
            break;
        }
      }
    }
    EndMode2D();
    EndDrawing();
  }

  // Don't re-run (probably never reached).
  return false;
}

int do_main() {
  InitWindow(600, 600, "Q");
  SetTargetFPS(60);
  while (show())
    ;
  CloseWindow();
  return 0;
}

#if defined(_WIN32)
int WinMain(void **_0, void **_1, void **_2, int _3) { return do_main(); }
#else
int main() { return do_main(); }
#endif
