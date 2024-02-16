#include <algorithm>
#include <array>
#include <deque>
#include <morton.h>
#include <random>
#include <raylib.h>
#include <vector>

std::array<float, 3> hsl2rgb(std::array<float, 3> hsl) {
  auto l = hsl[2], a = hsl[1] * std::min(l, 1.0f - l);
  auto min3 = [](float a, float b, float c) {
    return std::min(a, std::min(b, c));
  };
  auto f = [l, a, &hsl, &min3](float n) {
    auto k = std::remainder(n + hsl[0] * 12.0f, 12.0f);
    return l - a * std::max(-1.0f, min3(k - 3.0f, 9.0f - k, 1.0f));
  };
  return {f(0.0f), f(8.0f), f(4.0f)};
}

int do_main() {
  long constexpr SEED = 5432, N_PARTICLES = 10000,
                 MAX_N_QUEUE = N_PARTICLES / 4;
  float constexpr RADIUS = 0.00390625f, SATURATION = 1.0f, LIGHTNESS = 0.66f;
  uint64_t constexpr MASK = 0xffff'0000'0000'0000;

  std::mt19937 rng(SEED);
  std::uniform_real_distribution<float> udist;

  // Generate particles and then sort them by Morton order.
  std::vector<std::complex<float>> pp;
  for (long i = 0; i < N_PARTICLES; i++) {
    pp.emplace_back(udist(rng), udist(rng));
  }
  std::sort(pp.begin(), pp.end(), [](auto const &p, auto const &q) {
    return dyn::morton32(p) < dyn::morton32(q);
  });

  InitWindow(600, 600, "Q");
  SetTargetFPS(60);
  Camera2D cam{};
  cam.zoom = 600.0f;
  long q{};
  uint64_t prefix{};
  std::deque<Color> colors;
  while (!WindowShouldClose()) {
    BeginDrawing();
    BeginMode2D(cam);
    {
      ClearBackground(BLACK);
      // Color and then draw the particles.
      auto particle = pp[q];
      auto jpfx = dyn::morton32(particle) & MASK;
      if (jpfx != prefix) {
        auto hue = udist(rng);
        auto rgb = hsl2rgb({hue, SATURATION, LIGHTNESS});
        colors.emplace_front((unsigned char)(255.0f * rgb[0]),
                             (unsigned char)(255.0f * rgb[1]),
                             (unsigned char)(255.0f * rgb[2]), 255);
        if (colors.size() > MAX_N_QUEUE)
          colors.pop_back();
      } else {
        colors.push_front(colors.front());
        if (colors.size() > MAX_N_QUEUE)
          colors.pop_back();
      }
      prefix = jpfx;
      for (long n = 0; n < MAX_N_QUEUE; n++) {
        auto i = q - n;
        auto j = i < 0 ? i + N_PARTICLES : i;
        auto p = pp[j];
        auto c = n < colors.size() ? colors[n] : BLACK;
        c = Fade(c, 1.0f / (1.0f + std::pow(float(n), 0.25f)));
        DrawCircleV({p.real(), p.imag()}, RADIUS, c);
      }
      q = (q + 1) % N_PARTICLES;
    }
    EndMode2D();
    EndDrawing();
  }
  CloseWindow();

  return 0;
}

#if defined(_WIN32)
int WinMain(void **_0, void **_1, void **_2, int _3) { return do_main(); }
#else
int main() { return do_main(); }
#endif
