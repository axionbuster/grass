#include <algorithm>
#include <array>
#include <morton.h>
#include <random>
#include <raylib.h>
#include <vector>

int do_main() {
  long constexpr SEED = 1234, N_PARTICLES = 100, DIVISOR = 25,
                 C_QUOTIENT = N_PARTICLES / DIVISOR + !!(N_PARTICLES % DIVISOR);
  float constexpr RADIUS = 0.00390625f;

  std::array<Color, C_QUOTIENT> colors{};
  std::mt19937 rng(SEED);
  std::uniform_real_distribution<float> udist;
  // Generate a fixed palette of colors.
  for (auto &&c : colors) {
    struct {
      float r, g, b;
      [[nodiscard]] constexpr Color color() const {
        // r, g, b, a.
        return {(unsigned char)(255.0f * r), (unsigned char)(255.0f * g),
                (unsigned char)(255.0f * b), 255};
      }
    } a{udist(rng), udist(rng), udist(rng)};
    c = a.color();
  }

  // Generate particles and then sort them by Morton order.
  std::vector<std::complex<float>> pp;
  for (long i = 0; i < N_PARTICLES; i++) {
    pp.emplace_back(udist(rng), udist(rng));
  }
  std::sort(pp.begin(), pp.end(), [](auto const &p, auto const &q) {
    return dyn::morton32(p) < dyn::morton32(q);
  });

  struct Box {
    std::complex<float> ll, gg;
    [[nodiscard]] constexpr Rectangle rectangle() const {
      auto diff = gg - ll;
      return {ll.real(), ll.imag(), diff.real(), diff.imag()};
    }
  };
  std::vector<Box> boxes;
  for (long q = 0; q < C_QUOTIENT; q++) {
    auto ill = q * DIVISOR, igg = std::min((q + 1) * DIVISOR, N_PARTICLES) - 1L;
    boxes.emplace_back(pp[ill], pp[igg]);
  }

  InitWindow(600, 600, "Q");
  SetTargetFPS(90);
  Camera2D cam{};
  cam.zoom = 600.0f;
  while (!WindowShouldClose()) {
    BeginDrawing();
    BeginMode2D(cam);
    {
      ClearBackground(BLACK);
      // Color and then draw the particles.
      for (long i = 0; i < N_PARTICLES; i++) {
        auto color = colors[i / DIVISOR];
        auto particle = pp[i];
        DrawCircleV({particle.real(), particle.imag()}, RADIUS, color);
      }
      // Draw the bounding boxes.
      for (long i = 0; i < boxes.size(); i++) {
        auto color = colors[i];
        auto rectangle = boxes[i].rectangle();
        if (rectangle.width < 0) {
          rectangle.x += rectangle.width;
          rectangle.width = std::abs(rectangle.width);
        }
        if (rectangle.height < 0) {
          rectangle.y += rectangle.height;
          rectangle.height = std::abs(rectangle.height);
        }
        DrawRectangleLinesEx(rectangle, RADIUS, color);
      }
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
